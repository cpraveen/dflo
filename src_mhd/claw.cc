#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "claw.h"
#include "ic.h"

// Coefficients for SSP-RK scheme
double ark[3];
unsigned int n_rk;

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
                                       //const unsigned int degree,
                                       const FE_DGQArbitraryNodes<dim> &fe_scalar)
   :
   mpi_communicator (MPI_COMM_WORLD),
   triangulation(mpi_communicator),
   fe (fe_scalar, EulerEquations<dim>::n_components),
   dof_handler (triangulation),
   fe_cell (FE_DGQ<dim>(0)),
   dh_cell (triangulation),
   pcout (std::cout,(Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
   computing_timer (pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
{
   read_parameters (input_filename);
}

//------------------------------------------------------------------------------
// Constructor for Pk basis
//------------------------------------------------------------------------------
template <int dim>
ConservationLaw<dim>::ConservationLaw (const char *input_filename,
                                       //const unsigned int degree,
                                       const FE_DGP<dim> &fe_scalar)
:
   mpi_communicator (MPI_COMM_WORLD),
   triangulation(mpi_communicator),
   fe (fe_scalar, EulerEquations<dim>::n_components),
   dof_handler (triangulation),
   fe_cell (FE_DGQ<dim>(0)),
   dh_cell (triangulation),
   pcout (std::cout,(Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
   computing_timer (pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
{
   read_parameters (input_filename);
}

//------------------------------------------------------------------------------
// Read parameters from file and set some other parameters.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::read_parameters (const char *input_filename)
{
   AssertThrow(dim == 2, ExcMessage("Only works for 2-D"));
   
   ParameterHandler prm;
   Parameters::AllParameters<dim>::declare_parameters (prm);
   
   prm.read_input (input_filename);
   parameters.parse_parameters (prm);
   
   //pcout.set_condition (parameters.output == Parameters::Solver::verbose);
   
   // Create directory to save solution files
   if(Utilities::MPI::this_mpi_process(mpi_communicator)==0)
   {
      //system("mkdir -p output");
      
      // Save all parameters in xml format
      //std::ofstream xml_file ("input.xml");
      //prm.print_parameters (xml_file,  ParameterHandler::XML);
      
      // Save all parameters in txt format
      std::ofstream txt_file ("output/input.txt");
      prm.print_parameters (txt_file,  ParameterHandler::Text);
   }
   
   // Set coefficients for SSPRK
   if(fe.degree == 0)
   {
      ark[0] = 0.0;
      n_rk   = 1;
   }
   else if(fe.degree == 1)
   {
      ark[0] = 0.0;
      ark[1] = 1.0/2.0;
      n_rk   = 2;
   }
   else
   {
      ark[0] = 0.0;
      ark[1] = 3.0/4.0;
      ark[2] = 1.0/3.0;
      n_rk   = 3;
   }
}

//------------------------------------------------------------------------------
// Return mapping type based on selected type
//------------------------------------------------------------------------------
template <int dim>
const Mapping<dim,dim>& ConservationLaw<dim>::mapping() const
{
   if(parameters.mapping_type == Parameters::AllParameters<dim>::q1)
   {
      static MappingQ1<dim> m;
      return m;
   }
   else if(parameters.mapping_type == Parameters::AllParameters<dim>::q2)
   {
      static MappingQ<dim> m(2);
      return m;
   }
   else if(parameters.mapping_type == Parameters::AllParameters<dim>::cartesian)
   {
      static MappingCartesian<dim> m;
      return m;
   }
   else
   {
      AssertThrow (false, ExcMessage("Requested mapping type is unknown"));
      static MappingQ1<dim> m;
      return m;
   }

}

//------------------------------------------------------------------------------
// Compute dx, dy, dz for cartesian meshes.
// At present it only checks that dx == dy
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::compute_cartesian_mesh_size ()
{
   const double geom_tol = 1.0e-12;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      double xmin = 1.0e20, xmax = -1.0e20;
      double ymin = 1.0e20, ymax = -1.0e20;
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
         const Point<dim>& p = cell->face(f)->center();
         xmin = std::min(xmin, p[0]);
         xmax = std::max(xmax, p[0]);
         ymin = std::min(ymin, p[1]);
         ymax = std::max(ymax, p[1]);
      }
      double dx = xmax - xmin;
      double dy = ymax - ymin;
      AssertThrow(std::fabs(dx-dy) < geom_tol, ExcMessage("Cell is not square"));
   }
}

//------------------------------------------------------------------------------
// Our nodal basis uses Gauss points and we use Gauss quadrature with degree+1
// nodes. Then the mass matrix is diagonal. The mass matrix is integrated exactly
// for Cartesian cells but not exact for general cells.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::compute_inv_mass_matrix ()
{
   pcout << "Creating mass matrix ...\n";

   QGauss<dim> quadrature(fe.degree+1);
   unsigned int n_q_points = quadrature.size();
   FEValues<dim> fe_values(mapping(), fe, quadrature, update_values | update_JxW_values);
   
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   inv_mass_matrix.resize(triangulation.n_active_cells(),
                          Vector<double>(fe.dofs_per_cell));
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      unsigned int c = cell_number(cell);
      cell->get_dof_indices (dof_indices);
      fe_values.reinit(cell);
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
         inv_mass_matrix[c][i] = 0.0;
         for(unsigned int q=0; q<n_q_points; ++q)
            inv_mass_matrix[c][i] += fe_values.shape_value(i,q) *
                                     fe_values.shape_value(i,q) *
                                     fe_values.JxW(q);
         inv_mass_matrix[c][i] = 1.0 / inv_mass_matrix[c][i];
      }
   }
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
   TimerOutput::Scope t(computing_timer, "Setup");

   //pcout << "Allocating memory ...\n";
   
   dof_handler.clear();
   dof_handler.distribute_dofs (fe);
   
   locally_owned_dofs = dof_handler.locally_owned_dofs ();
   DoFTools::extract_locally_relevant_dofs (dof_handler,
                                            locally_relevant_dofs);
   
   // Size all of the fields.
   current_solution.reinit  (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
   right_hand_side.reinit   (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
   old_solution.reinit 		(locally_owned_dofs, mpi_communicator);
   predictor.reinit 		(locally_owned_dofs, mpi_communicator);
   newton_update.reinit 	(locally_owned_dofs, mpi_communicator);
   
   cell_average.resize 		(triangulation.n_active_cells(),
							       Vector<double>(EulerEquations<dim>::n_components));
   
   mu_shock.reinit 			(triangulation.n_active_cells());
   shock_indicator.reinit 	(triangulation.n_active_cells());
   jump_indicator.reinit 	(triangulation.n_active_cells());
   
   pcout << std::endl
         << "   Number of active cells:       " << triangulation.n_global_active_cells()
         << std::endl
         << "   Number of degrees of freedom: " << dof_handler.n_dofs()
         << std::endl << std::endl;

   // create map from (level,index) to cell number
   unsigned int index=0;
   
   for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
		cell!=triangulation.end(); ++cell, ++index)
			cell->set_user_index(index);

   // compute mass matrix
   compute_inv_mass_matrix ();

   if(parameters.mapping_type == Parameters::AllParameters<dim>::cartesian)
      compute_cartesian_mesh_size ();

   if(parameters.mapping_type != Parameters::AllParameters<dim>::cartesian)
      return;

   // For each cell, find neighbourig cell
   // This is needed for limiter
   // CHECK: Should the size be n_active_cells() ?
   lcell.resize(triangulation.n_active_cells());
   rcell.resize(triangulation.n_active_cells());
   bcell.resize(triangulation.n_active_cells());
   tcell.resize(triangulation.n_active_cells());

   //const double EPS = 1.0e-10;
   typename DoFHandler<dim>::active_cell_iterator
      cell = dh_cell.begin_active(),
      endc = dh_cell.end();
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      unsigned int c = cell_number(cell);
      lcell[c] = endc;
      rcell[c] = endc;
      bcell[c] = endc;
      tcell[c] = endc;
      double dx = cell->diameter() / std::sqrt(1.0*dim);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
         if (! cell->at_boundary(face_no))
         {
            const typename DoFHandler<dim>::cell_iterator
               neighbor = cell->neighbor(face_no);
            Assert(neighbor->level() == cell->level() || neighbor->level() == cell->level()-1,
                   ExcInternalError());
            Tensor<1,dim, double> dr = neighbor->center() - cell->center();
	    //Point<dim> dr = neighbor->center() - cell->center();
            if(dr[0] < -0.5*dx)
               lcell[c] = neighbor;
            else if(dr[0] > 0.5*dx)
               rcell[c] = neighbor;
            else if(dr[1] < -0.5*dx)
               bcell[c] = neighbor;
            else if(dr[1] > 0.5*dx)
               tcell[c] = neighbor;
            else
            {
               std::cout << "Did not find all neighbours\n";
               std::cout << "dx, dy = " << dr[0] << "  " << dr[1] << std::endl;
               //std::cout << "dx, dy = " << dr(0) << "  " << dr(1) << std::endl;
               exit(0);
            }
         }
   }
}

//------------------------------------------------------------------------------
// Create mesh worker for explicit integration
// This computes only the right hand side
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::setup_mesh_worker (IntegratorExplicit<dim>& integrator)
{
   pcout << "Setting up mesh worker ...\n";

   const unsigned int n_gauss_points = fe.degree + 1;
   integrator.info_box.initialize_gauss_quadrature(n_gauss_points,
                                                   n_gauss_points,
                                                   n_gauss_points);
   
   integrator.info_box.initialize_update_flags   ();
   integrator.info_box.add_update_flags_all 	    (update_values | update_JxW_values);
   integrator.info_box.add_update_flags_cell     (update_gradients);
   integrator.info_box.add_update_flags_boundary (update_normal_vectors | update_quadrature_points);
   integrator.info_box.add_update_flags_face     (update_normal_vectors);
   
   integrator.info_box.initialize (fe, mapping());
   
   //NamedData< LA::Vector<double>* > rhs;
   AnyData rhs;
   LA::Vector<double>* data = &right_hand_side;
   rhs.add< LA::Vector<double>* > (data, "RHS");
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
   TimerOutput::Scope t(computing_timer, "Time step");

   // No need to compute time step for stationary flows
   if(parameters.is_stationary == true)
      return;

   // Update size of dt array since adaptation might have been performed
   dt.reinit (triangulation.n_active_cells());

   // If time step given in input file, then use it. This is only for global time stepping
   if(parameters.time_step_type == "global" && parameters.cfl <= 0.0)
   {
      dt = parameters.time_step;
      return;
   }

   if(parameters.mapping_type == Parameters::AllParameters<dim>::cartesian)
      compute_time_step_cartesian();
   else
      compute_time_step_q();


   // For global time step, use the minimum value
   if(parameters.time_step_type == "global")
   {
      if(global_dt > 0 && parameters.time_step > 0)
         global_dt = std::min(global_dt, parameters.time_step);
      if(elapsed_time + global_dt > parameters.final_time)
         global_dt = parameters.final_time - elapsed_time;
      dt = global_dt;
   }

}

//------------------------------------------------------------------------------
// Compute local time step for each cell
// This function is specified to cartesian cells.
//------------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::compute_time_step_cartesian ()
{
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   global_dt = 1.0e20;
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      const unsigned int c = cell_number (cell);
      const double h = cell->diameter() / std::sqrt(1.0*dim);
      const double sonic = EulerEquations<dim>::sound_speed(cell_average[c]);
      const double density = cell_average[c][EulerEquations<dim>::density_component];
      
      double max_eigenvalue = 0.0;
      for (unsigned int d=0; d<dim; ++d)
         max_eigenvalue += (sonic + std::fabs(cell_average[c][d]/density))/h;
      
      dt(c) = parameters.cfl / max_eigenvalue / (2.0*fe.degree + 1.0);
      
      global_dt = std::min(global_dt, dt(c));
   }
   
   global_dt = -global_dt;
   global_dt = -Utilities::MPI::max (global_dt, mpi_communicator);
}

//------------------------------------------------------------------------------
// Compute local time step for each cell
// We compute speed at quadrature points and take maximum over these values.
// This speed is used to compute time step in each cell.
//------------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::compute_time_step_q ()
{
   QIterated<dim>   quadrature_formula(QTrapez<1>(), 3);
   const unsigned int   n_q_points = quadrature_formula.size();
   
   FEValues<dim> fe_values (mapping(), fe,
                            quadrature_formula,
                            update_values);
   std::vector<Vector<double> > solution_values(n_q_points,
                                                Vector<double>(EulerEquations<dim>::n_components));

                                                   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   global_dt = 1.0e20;
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
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
      dt(c) = parameters.cfl * h / max_eigenvalue / (2.0*fe.degree + 1.0) / dim;
      
      global_dt = std::min(global_dt, dt(c));
   }
   global_dt = -global_dt;
   global_dt = -Utilities::MPI::max (global_dt, mpi_communicator);
}

//------------------------------------------------------------------------------
// Compute cell average solution
//------------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::compute_cell_average ()
{
   TimerOutput::Scope t(computing_timer, "Cell average");

   QGauss<dim>   quadrature_formula(fe.degree+1);
   const unsigned int n_q_points = quadrature_formula.size();
   
   FEValues<dim> fe_values (mapping(), fe,
                            quadrature_formula,
                            update_values | update_JxW_values);
   std::vector<Vector<double> > solution_values(n_q_points,
                                                Vector<double>(EulerEquations<dim>::n_components));
                                                
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   // compute cell average for ghost cells also.
   for (; cell!=endc; ++cell)
   if(!cell->is_artificial())
   {
      unsigned int cell_no = cell_number(cell);
      {
         fe_values.reinit (cell);
         fe_values.get_function_values (current_solution, solution_values);
         
         cell_average[cell_no] = 0.0;
         
         for (unsigned int q=0; q<n_q_points; ++q)
            for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
               cell_average[cell_no][c] += solution_values[q][c] * fe_values.JxW(q);
         
         cell_average[cell_no] /= cell->measure();
      }
   }
}

//------------------------------------------------------------------------------
// Computes total angular momentum in the computational domain
//------------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::compute_angular_momentum ()
{
   AssertThrow(dim==2, ExcNotImplemented());
   
   QGauss<dim>   quadrature_formula(fe.degree+1);
   const unsigned int n_q_points = quadrature_formula.size();
   
   FEValues<dim> fe_values (mapping(), fe,
                            quadrature_formula,
                            update_values | update_q_points | update_JxW_values);
   const FEValuesExtractors::Vector momentum (0);
   std::vector< Tensor<1,dim> > momentum_values(n_q_points);
   
   double angular_momentum = 0;
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      fe_values.reinit(cell);
      fe_values[momentum].get_function_values(current_solution, momentum_values);
      const std::vector<Point<dim> >& qp = fe_values.get_quadrature_points();
      for(unsigned int q=0; q<n_q_points; ++q)
      {
         double cross = qp[q][0] * momentum_values[q][1] - qp[q][1] * momentum_values[q][0];
         angular_momentum += cross * fe_values.JxW(q);
      }
   }
   
   angular_momentum = Utilities::MPI::sum (angular_momentum, mpi_communicator);
   pcout << "Total angular momentum: "
         << elapsed_time << "  " <<  angular_momentum << std::endl;
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
ConservationLaw<dim>::solve (LA::Vector<double> &newton_update)// , double              current_residual)
{
   TimerOutput::Scope t(computing_timer, "Solve");
   
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
      if(cell->is_locally_owned())
      {
         const unsigned int cell_no = cell_number (cell);
         
         cell->get_dof_indices (dof_indices);
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            newton_update(dof_indices[i]) = dt(cell_no) *
                                            right_hand_side(dof_indices[i]) *
                                            inv_mass_matrix[cell_no][i];
      }
   return std::pair<unsigned int, double> (0,0);
}

//------------------------------------------------------------------------------
// Perform RK update
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::iterate_explicit (IntegratorExplicit<dim>& integrator,
                                             LA::Vector<double>& newton_update,
                                             double& res_norm0, double& res_norm)
{
   
   // Loop for newton iterations or RK stages
   for(unsigned int rk=0; rk<n_rk; ++rk)
   {
      // set time in boundary condition
      // NOTE: We need to check if this is time accurate.
      for (unsigned int boundary_id=0; boundary_id<Parameters::AllParameters<dim>::max_n_boundaries;
           ++boundary_id)
      {
         parameters.boundary_conditions[boundary_id].values.set_time(elapsed_time);
      }
      
      assemble_system (integrator);
      
      res_norm = right_hand_side.l2_norm();
      if(rk == 0) res_norm0 = res_norm;
      
      //std::pair<unsigned int, double> convergence = solve (newton_update, res_norm);
      
      {
         TimerOutput::Scope t(computing_timer, "RK update");
         
         // current_solution = (1-ark)*current_solution + ark*old_solution + (1-ark)*newton_update
         current_solution.sadd(1.0-ark[rk], ark[rk], old_solution, 1-ark[rk], newton_update);
      }
      
      compute_cell_average ();
      compute_shock_indicator ();
      apply_limiter ();
   }
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
   
//   {
//      // Use this refine around corner in forward step
//      unsigned int nrefine = 2;
//      for(unsigned int i=0; i<nrefine; ++i)
//         refine_forward_step ();
//   }
   
   /*
   static const HyperBallBoundary<dim> boundary_description;
   triangulation.set_boundary (1, boundary_description);
   */
   
   setup_system();
   set_initial_condition ();
   
   // Refine the initial mesh
   if (parameters.do_refine == true)
      for (unsigned int i=0; i<parameters.shock_levels; ++i)
      {
         Vector<double> refinement_indicators (triangulation.n_active_cells());
         compute_refinement_indicators(refinement_indicators);
         refine_grid(refinement_indicators);
         set_initial_condition ();
      }
   
   // Cell average of initial condition
   compute_cell_average ();

   // Limit the initial condition
   compute_shock_indicator ();
   apply_limiter();
   old_solution = current_solution;
   
   // Reset time/iteration counters
   elapsed_time = 0;
   time_iter = 0;
   
   // Save initial condition to file
   output_results ();
   
   // We then enter into the main time stepping loop.

   // Setup variables to decide if we need to save solution 
   double next_output_time = elapsed_time + parameters.output_time_step;
   int next_output_iter = time_iter + parameters.output_iter_step;

   // Variable to control grid refinement
   double next_refine_time = elapsed_time + parameters.refine_time_step;
   int    next_refine_iter = time_iter + parameters.refine_iter_step;

   std::vector<double> residual_history;
   IntegratorExplicit<dim> integrator_explicit (dof_handler);
   setup_mesh_worker (integrator_explicit);
   
   while (elapsed_time < parameters.final_time)
   {
      // compute time step in each cell using cfl condition
      compute_time_step ();
      
      pcout << std::scientific << "It = " << std::setw(6) << time_iter+1
            << ", T = " << std::setprecision(4) << elapsed_time + global_dt
            << ", dt = " << std::setprecision(4) << global_dt
            << ", cfl = " << std::setprecision(2) << parameters.cfl
            << std::endl;
      
      //unsigned int nonlin_iter = 0;
      double res_norm0 = 1.0;
      double res_norm  = 1.0;

      iterate_explicit(integrator_explicit, newton_update, res_norm0, res_norm);
      
      // Update counters
      elapsed_time += global_dt;
      ++time_iter;
      
      if(time_iter % parameters.ang_mom_step == 0)
         compute_angular_momentum();

      residual_history.push_back (res_norm);
      
      // Save solution for visualization
      if (elapsed_time >= next_output_time || time_iter == next_output_iter 
            || std::fabs(elapsed_time-parameters.final_time) < 1.0e-13)
      {
         output_results ();
         next_output_time = elapsed_time + parameters.output_time_step;
         next_output_iter = time_iter + parameters.output_iter_step;
      }
      
      old_solution = current_solution;
      
      if (parameters.do_refine == true && 
          (elapsed_time >= next_refine_time || time_iter == next_refine_iter))
      {
         Vector<double> refinement_indicators (triangulation.n_active_cells());
         compute_refinement_indicators(refinement_indicators);
         
         refine_grid(refinement_indicators);
         
         //newton_update.reinit (locally_owned_dofs, mpi_communicator);

         next_refine_time = elapsed_time + parameters.refine_time_step;
         next_refine_iter = time_iter + parameters.refine_iter_step;

         // We may need to reduce the cfl after refinement, only for steady state
         // problems when using gmres
         //parameters.cfl = 1.2;
      }
   }
   
   computing_timer.print_summary ();
   computing_timer.reset ();
   pcout << std::endl;
}

template class ConservationLaw<2>;
