#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>

#include <lac/vector.h>
#include <lac/compressed_sparsity_pattern.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_in.h>
#include <grid/tria_boundary_lib.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>

#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>

#include <numerics/vector_tools.h>
#include <numerics/solution_transfer.h>
#include <numerics/matrix_tools.h>

#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_precondition.h>
#include <lac/trilinos_solver.h>

#include <lac/solver_gmres.h>
#include <lac/sparse_direct.h>
#include <lac/precondition_block.h>

#include <Sacado.hpp>


#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "claw.h"

// Coefficients for 2-stage RK scheme
const double ark[3] = {0.0, 3.0/4.0, 1.0/3.0};

using namespace dealii;

//------------------------------------------------------------------------------
// @sect4{ConservationLaw::ConservationLaw}
//
// There is nothing much to say about
// the constructor. Essentially, it
// reads the input file and fills the
// parameter object with the parsed
// values:
//------------------------------------------------------------------------------
template <int dim>
ConservationLaw<dim>::ConservationLaw (const char *input_filename,
                                       const unsigned int mapping,
                                       const unsigned int degree)
   :
   mapping (mapping),
   fe (FE_DGQ<dim>(degree), EulerEquations<dim>::n_components),
   dof_handler (triangulation),
   fe0 (FE_DGQ<dim>(0), EulerEquations<dim>::n_components),
   dof_handler0 (triangulation),
   fe_cell (FE_DGQ<dim>(0)),
   dh_cell (triangulation),
   verbose_cout (std::cout, false)
{
   ParameterHandler prm;
   Parameters::AllParameters<dim>::declare_parameters (prm);
   
   prm.read_input (input_filename);
   parameters.parse_parameters (prm);
   
   verbose_cout.set_condition (parameters.output == Parameters::Solver::verbose);

   // Save all parameters in xml format
   std::ofstream xml_file ("input.xml");
   prm.print_parameters (xml_file,  ParameterHandler::XML);
}

//------------------------------------------------------------------------------
// @sect4{ConservationLaw::setup_system}
//
// The following (easy) function is called
// each time the mesh is changed. All it
// does is to resize the Trilinos matrix
// according to a sparsity pattern that we
// generate as in all the previous tutorial
// programs.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::setup_system ()
{
   //DoFRenumbering::Cuthill_McKee (dof_handler);
   
   dof_handler.clear();
   dof_handler.distribute_dofs (fe);
   
   // Size all of the fields.
   old_solution.reinit (dof_handler.n_dofs());
   current_solution.reinit (dof_handler.n_dofs());
   predictor.reinit (dof_handler.n_dofs());
   right_hand_side.reinit (dof_handler.n_dofs());
   
   dof_handler0.clear();
   dof_handler0.distribute_dofs (fe0);
   cell_average.reinit (dof_handler0.n_dofs());
   
   // Used for cell data like time step
   dh_cell.clear();
   dh_cell.distribute_dofs (fe_cell);
   mu_shock.reinit (dh_cell.n_dofs());

   if(parameters.solver == Parameters::Solver::rk3)
   {
      std::cout << "Creating mass matrix ...\n";
      CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
      sparsity_pattern.copy_from(c_sparsity);
      system_matrix.reinit (sparsity_pattern);
      MatrixCreator::create_mass_matrix (mapping, 
                                         dof_handler, 
                                         QGauss<dim>(fe.degree+1), 
                                         system_matrix);
      // Compute the inverse of the blocks
      inv_mass_matrix.initialize (system_matrix, fe.dofs_per_cell);
      // Once inv_mass_matrix is created, we can delete system_matrix
      system_matrix.clear();
   }
   else
   {
      CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
      DoFTools::make_flux_sparsity_pattern (dof_handler, c_sparsity);
      sparsity_pattern.copy_from(c_sparsity);

      system_matrix.reinit (sparsity_pattern);
   }

   // For each cell, find neighbourig cell
   // This is needed for limiter
   lcell.resize(triangulation.n_cells());
   rcell.resize(triangulation.n_cells());
   bcell.resize(triangulation.n_cells());
   tcell.resize(triangulation.n_cells());

   const double EPS = 1.0e-10;
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler0.begin_active(),
      endc = dof_handler0.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      lcell[c] = endc;
      rcell[c] = endc;
      bcell[c] = endc;
      tcell[c] = endc;
      double dx = cell->diameter() / std::sqrt(2.0);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
         if (! cell->at_boundary(face_no))
         {
            const typename DoFHandler<dim>::cell_iterator
               neighbor = cell->neighbor(face_no);
            Assert(neighbor->level() == cell->level() || neighbor->level() == cell->level()-1,
                   ExcInternalError());
            Point<dim> dr = neighbor->center() - cell->center();
            if(dr(0) < -0.5*dx)
               lcell[c] = neighbor;
            else if(dr(0) > 0.5*dx)
               rcell[c] = neighbor;
            else if(dr(1) < -0.5*dx)
               bcell[c] = neighbor;
            else if(dr(1) > 0.5*dx)
               tcell[c] = neighbor;
            else
            {
               std::cout << "Did not find all neighbours\n";
               std::cout << "dx, dy = " << dr(0) << "  " << dr(1) << std::endl;
               exit(0);
            }
         }
   }

   // Visualize sparsity pattern
   //std::ofstream out ("sparsity_pattern.1");
   //sparsity_pattern.print_gnuplot (out);
   //abort ();
}

//------------------------------------------------------------------------------
// Create mesh worker for implicit integration
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::setup_mesh_worker (IntegratorImplicit<dim>& integrator)
{   
   const unsigned int n_gauss_points = fe.degree + 1;
   integrator.info_box.initialize_gauss_quadrature(n_gauss_points,
                                                   n_gauss_points,
                                                   n_gauss_points);
   
   integrator.info_box.initialize_update_flags ();
   integrator.info_box.add_update_flags_all (update_values | 
                                             update_quadrature_points |
                                             update_JxW_values);
   integrator.info_box.add_update_flags_cell     (update_gradients);
   integrator.info_box.add_update_flags_boundary (update_normal_vectors | update_gradients);
   integrator.info_box.add_update_flags_face     (update_normal_vectors | update_gradients);
   
   integrator.info_box.initialize (fe, mapping);
   
   integrator.assembler.initialize (system_matrix, right_hand_side);
}

//------------------------------------------------------------------------------
// Create mesh worker for explicit integration
// This computes only the right hand side
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::setup_mesh_worker (IntegratorExplicit<dim>& integrator)
{   
   const unsigned int n_gauss_points = fe.degree + 1;
   integrator.info_box.initialize_gauss_quadrature(n_gauss_points,
                                                   n_gauss_points,
                                                   n_gauss_points);
   
   integrator.info_box.initialize_update_flags ();
   integrator.info_box.add_update_flags_all (update_values | 
                                             update_quadrature_points |
                                             update_JxW_values);
   integrator.info_box.add_update_flags_cell     (update_gradients);
   integrator.info_box.add_update_flags_boundary (update_normal_vectors | update_gradients); // TODO:ADIFF
   integrator.info_box.add_update_flags_face     (update_normal_vectors | update_gradients); // TODO:ADIFF
   
   integrator.info_box.initialize (fe, mapping);
   
   NamedData< Vector<double>* > rhs;
   Vector<double>* data = &right_hand_side;
   rhs.add (data, "RHS");
   integrator.assembler.initialize (rhs);
}

//------------------------------------------------------------------------------
// Compute local time step for each cell
// We compute speed at quadrature points and take maximum over these values.
// This speed is used to compute time step in each cell.
//------------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::compute_time_step ()
{
   // No need to compute time step for stationary flows
   if(parameters.is_stationary == true)
      return;

   // Update size of dt array since adaptation might have been performed
   dt.reinit (triangulation.n_cells());

   // If time step given in input file, then use it. This is only for global time stepping
   if(parameters.time_step_type == "global" && parameters.time_step > 0.0)
   {
      dt = parameters.time_step;
      return;
   }


   //QGaussLobatto<dim>   quadrature_formula(fe.degree+1);
   QIterated<dim>   quadrature_formula(QTrapez<1>(), 3);
   const unsigned int   n_q_points = quadrature_formula.size();
   
   FEValues<dim> fe_values (fe, 
                            quadrature_formula,
                            update_values);
   std::vector<Vector<double> > solution_values(n_q_points,
                                                Vector<double>(EulerEquations<dim>::n_components));
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   double min_time_step = 1.0e20;

   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      fe_values.get_function_values (current_solution, solution_values);
      
      double max_eigenvalue = 0.0;
      for (unsigned int q=0; q<n_q_points; ++q)
      {
         max_eigenvalue = std::max(max_eigenvalue,
                                   EulerEquations<dim>::max_eigenvalue (solution_values[q]));
      }
      
      const unsigned int c = cell_number (cell);
      const double h = cell->diameter() / std::sqrt(1.0*dim);
      dt(c) = parameters.cfl * h / max_eigenvalue / (2.0*fe.degree + 1.0);

      min_time_step = std::min(min_time_step, dt(c));
   }


   // For global time step, use the minimum value
   if(parameters.time_step_type == "global")
   {
      dt = min_time_step;
      parameters.time_step = min_time_step;
   }

}

//------------------------------------------------------------------------------
// Compute cell average solution
//------------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::compute_cell_average ()
{
   QGauss<dim>   quadrature_formula(fe.degree+1);
   const unsigned int n_q_points = quadrature_formula.size();
   
   FEValues<dim> fe_values (fe,
                            quadrature_formula,
                            update_values | update_JxW_values);
   std::vector<Vector<double> > solution_values(n_q_points,
                                                Vector<double>(EulerEquations<dim>::n_components));
   std::vector<unsigned int> cell_dofs(fe0.dofs_per_cell);
   
   cell_average = 0;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell0= dof_handler0.begin_active();
   
   for (; cell!=endc; ++cell, ++cell0)
   {
      fe_values.reinit (cell);
      fe_values.get_function_values (current_solution, solution_values);
      
      cell0->get_dof_indices(cell_dofs);
      
      for (unsigned int q=0; q<n_q_points; ++q)
         for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
            cell_average(cell_dofs[c]) += solution_values[q][c] * fe_values.JxW(q);
      
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         cell_average(cell_dofs[c]) /= cell->measure();
      
   }
   
}

//------------------------------------------------------------------------------
// @sect4{ConservationLaw::solve}
//
// Here, we actually solve the linear system,
// using either of Trilinos' Aztec or Amesos
// linear solvers. The result of the
// computation will be written into the
// argument vector passed to this
// function. The result is a pair of number
// of iterations and the final linear
// residual.
//------------------------------------------------------------------------------
template <int dim>
std::pair<unsigned int, double>
ConservationLaw<dim>::solve (Vector<double> &newton_update, 
                             double          current_residual)
{
   newton_update = 0;

   switch (parameters.solver)
   {
      case Parameters::Solver::umfpack:
      {
         SparseDirectUMFPACK  solver;
         solver.initialize(system_matrix);
         solver.vmult (newton_update, right_hand_side);
         return std::pair<unsigned int, double> (1, 0);
      }

      case Parameters::Solver::gmres:
      {
         SolverControl   solver_control (parameters.max_iterations, 
                                         parameters.linear_residual * current_residual);
         SolverGMRES<>::AdditionalData  gmres_data (30, true, true);
         SolverGMRES<>   solver (solver_control, gmres_data);

         PreconditionBlockSSOR<SparseMatrix<double>, float> preconditioner;
         preconditioner.initialize(system_matrix, fe.dofs_per_cell);

         solver.solve (system_matrix, 
                       newton_update, 
                       right_hand_side,
                       preconditioner);

         return std::pair<unsigned int, double> (solver_control.last_step(),
                                                 solver_control.last_value());
      }

      // We have equation M*du/dt = rhs, where M = mass matrix
      case Parameters::Solver::rk3:
      {
         inv_mass_matrix.vmult (newton_update, right_hand_side);

         // Multiply newton_update by time step dt
         std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
         typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
         for (; cell!=endc; ++cell)
         {
            const unsigned int cell_no = cell_number (cell);

            cell->get_dof_indices (dof_indices);
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
               newton_update(dof_indices[i]) *= dt(cell_no);
         }
         return std::pair<unsigned int, double> (0,0);
      }

      default:
         Assert (false, ExcNotImplemented());
   }
   
   return std::pair<unsigned int, double> (0,0);
}
//------------------------------------------------------------------------------
// @sect4{ConservationLaw::run}

// This function contains the top-level logic
// of this program: initialization, the time
// loop, and the inner Newton iteration.
//
// At the beginning, we read the mesh file
// specified by the parameter file, setup the
// DoFHandler and various vectors, and then
// interpolate the given initial conditions
// on this mesh. We then perform a number of
// mesh refinements, based on the initial
// conditions, to obtain a mesh that is
// already well adapted to the starting
// solution. At the end of this process, we
// output the initial solution.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::run ()
{
   {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      
      std::ifstream input_file(parameters.mesh_filename.c_str());
      Assert (input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));
      
      if(parameters.mesh_type == "ucd")
         grid_in.read_ucd(input_file);
      else if(parameters.mesh_type == "gmsh")
         grid_in.read_msh(input_file);
   }

   /*
   static const HyperBallBoundary<dim> boundary_description;
   triangulation.set_boundary (1, boundary_description);
   */
   
   setup_system();
   
   VectorTools::interpolate(dof_handler,
                            parameters.initial_conditions, old_solution);
   current_solution = old_solution;
   predictor = old_solution;
   
   // Refine the initial mesh
   if (parameters.do_refine == true)
      for (unsigned int i=0; i<parameters.shock_levels; ++i)
      {
         Vector<double> refinement_indicators (triangulation.n_active_cells());
         
         compute_refinement_indicators(refinement_indicators);
         refine_grid(refinement_indicators);
         
         VectorTools::interpolate(dof_handler,
                                  parameters.initial_conditions, old_solution);
         current_solution = old_solution;
         predictor = old_solution;
      }
   
   output_results ();
   
   // We then enter into the main time
   // stepping loop. 
   Vector<double> newton_update (dof_handler.n_dofs());
   
   double time = 0;
   int time_iter = 0;

   // Setup variables to decide if we need to save solution 
   double next_output_time = time + parameters.output_time_step;
   int next_output_iter = time_iter + parameters.output_iter_step;

   // Variable to control grid refinement
   double next_refine_time = time + parameters.refine_time_step;
   int    next_refine_iter = time_iter + parameters.refine_iter_step;

   std::vector<double> residual_history;
   
   while (time < parameters.final_time)
   {
      std::cout << std::endl;
      std::cout << "T=" << time << ", dt=" << parameters.time_step
                << ", cfl=" << parameters.cfl << std::endl
		          << "   Number of active cells:       "
		          << triangulation.n_active_cells()
		          << std::endl
		          << "   Number of degrees of freedom: "
		          << dof_handler.n_dofs()
		          << std::endl
		          << std::endl;
      
      std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
                << "   _____________________________________" << std::endl;
      
      // compute time step in each cell using cfl condition
      compute_time_step ();

      unsigned int nonlin_iter = 0;
      double res_norm0 = 1.0;
      double res_norm  = 1.0;
      
      // With global time stepping, we can use predictor as initial
      // guess for the implicit scheme.
      if(parameters.time_step_type == "global")
         current_solution = predictor;
      
      IntegratorImplicit<dim> integrator_implicit (dof_handler);
      setup_mesh_worker (integrator_implicit);
      IntegratorExplicit<dim> integrator_explicit (dof_handler);
      setup_mesh_worker (integrator_explicit);

      // Loop for newton iterations or RK stages
      while(true)
      {
         if(parameters.solver == Parameters::Solver::rk3)
            assemble_system (integrator_explicit);
         else
            assemble_system (integrator_implicit);
         
         res_norm = right_hand_side.l2_norm();
         if(nonlin_iter == 0) res_norm0 = res_norm;

         std::pair<unsigned int, double> convergence
            = solve (newton_update, res_norm);
            
         current_solution += newton_update;

         if(parameters.solver == Parameters::Solver::rk3)
         {
            // current_solution = ark*old_solution + (1-ark)*current_solution
            current_solution.sadd (1.0-ark[nonlin_iter], ark[nonlin_iter], old_solution);
         }
         
         compute_cell_average ();
         apply_limiter ();
            
         std::printf("   %-16.3e %04d        %-5.2e\n",
                     res_norm, convergence.first, convergence.second);
         
         ++nonlin_iter;

         // Check that newton iterations converged
         if(parameters.solver == Parameters::Solver::gmres &&
            nonlin_iter == parameters.max_nonlin_iter &&
            std::fabs(res_norm) > 1.0e-10)
            AssertThrow (nonlin_iter <= parameters.max_nonlin_iter, 
                         ExcMessage ("No convergence in nonlinear solver"));

         // check stopping criterion
         if(parameters.solver == Parameters::Solver::rk3 &&
            nonlin_iter == 3) // 3-stage RK
            break;
         else if(parameters.solver == Parameters::Solver::gmres &&
                 (nonlin_iter == parameters.max_nonlin_iter ||
                  std::fabs(res_norm) <= 1.0e-10))
            break;
         else if(parameters.solver == Parameters::Solver::umfpack)
            break;
      }
      
      // Update counters
      time += parameters.time_step;
      ++time_iter;

      // Increase cfl
      if(parameters.solver == Parameters::Solver::gmres && 
         parameters.time_step_type == "local" &&
         time_iter >= 2)
      {
         //parameters.cfl *= 1.2;
         double factor = residual_history.back() / res_norm;
         factor = std::min( factor, 2.0 );
         factor = std::max( factor, 0.5 );
         parameters.cfl *= factor;
      }

      residual_history.push_back (res_norm);
      
      // Save solution for visualization
      if (time >= next_output_time || time_iter == next_output_iter)
      {
         output_results ();
         next_output_time = time + parameters.output_time_step;
         next_output_iter = time_iter + parameters.output_iter_step;
      }
      
      // Compute predictor only for global time stepping
      // For local time stepping, this is meaningless
      if(parameters.time_step_type == "global")
      {
         predictor = current_solution;
         predictor.sadd (2.0, -1.0, old_solution);
      }
      
      old_solution = current_solution;
      
      if (parameters.do_refine == true && 
          (time >= next_refine_time || time_iter == next_refine_iter))
      {
         Vector<double> refinement_indicators (triangulation.n_active_cells());
         compute_refinement_indicators(refinement_indicators);
         
         refine_grid(refinement_indicators);
         
         newton_update.reinit (dof_handler.n_dofs());

         next_refine_time = time + parameters.refine_time_step;
         next_refine_iter = time_iter + parameters.refine_iter_step;

         // We may need to reduce the cfl after refinement, only for steady state
         // problems when using gmres
         //parameters.cfl = 1.2;
      }
   }
}

template class ConservationLaw<2>;
