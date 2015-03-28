#ifndef __CLAW_H__
#define __CLAW_H__

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>
#include <base/timer.h>

#include <lac/vector.h>
#include <lac/parallel_vector.h>

#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>

#include <distributed/tria.h>
#include <distributed/grid_refinement.h>
#include <distributed/solution_transfer.h>

#include <vector>

#include "parameters.h"
#include "integrator.h"

using namespace dealii;

namespace LA
{
   using namespace ::parallel::distributed;
}

//-----------------------------------------------------------------------------
// Data needed for positivity limiter
//-----------------------------------------------------------------------------
template <int dim>
struct PosLimData
{
   PosLimData(const dealii::FESystem<dim>    &fe,
              const dealii::Mapping<dim,dim> &mapping,
              const std::pair<unsigned int,unsigned int> &local_range);
   unsigned int ngll;
   
   dealii::Quadrature<dim> quadrature_x;
   dealii::Quadrature<dim> quadrature_y;

   unsigned int n_q_points;
   
   dealii::FEValues<dim> fe_values_x;
   dealii::FEValues<dim> fe_values_y;

   std::vector<double> density_values, energy_values;
   std::vector< Tensor<1,dim> > momentum_values;

   std::vector<unsigned int> local_dof_indices;
   std::pair<unsigned int, unsigned int> local_range;
};

template <int dim>
PosLimData<dim>::PosLimData(const dealii::FESystem<dim>    &fe,
                            const dealii::Mapping<dim,dim> &mapping,
                            const std::pair<unsigned int,unsigned int> &local_range)
:
   ngll ((fe.degree+3)%2==0 ? (fe.degree+3)/2 : (fe.degree+4)/2),
   quadrature_x (QGaussLobatto<1>(ngll), QGauss<1>(fe.degree+1)),
   quadrature_y (QGauss<1>(fe.degree+1), QGaussLobatto<1>(ngll)),
   n_q_points (quadrature_x.size()),
   fe_values_x (mapping, fe, quadrature_x, update_values),
   fe_values_y (mapping, fe, quadrature_y, update_values),
   density_values (n_q_points),
   energy_values (n_q_points),
   momentum_values (n_q_points),
   local_dof_indices (fe.dofs_per_cell),
   local_range (local_range)
{
   
}

// @sect3{Conservation law class}

// Here finally comes the class that
// actually does something with all
// the Euler equation and parameter
// specifics we've defined above. The
// public interface is pretty much
// the same as always (the
// constructor now takes the name of
// a file from which to read
// parameters, which is passed on the
// command line). The private
// function interface is also pretty
// similar to the usual arrangement,
// with the
// <code>assemble_system</code>
// function split into three parts:
// one that contains the main loop
// over all cells and that then calls
// the other two for integrals over
// cells and faces, respectively.
template <int dim>
class ConservationLaw
{
public:
   ConservationLaw (const char *input_filename,
                    const unsigned int degree,
                    const dealii::FE_DGQArbitraryNodes<dim> &fe_scalar);
   ConservationLaw (const char *input_filename,
                    const unsigned int degree,
                    const dealii::FE_DGP<dim> &fe_scalar);
   void run ();
   
private:

   void read_parameters (const char *file_name);
   const Mapping<dim,dim>& mapping() const;
   void compute_cartesian_mesh_size ();
   void compute_inv_mass_matrix();
   void setup_system ();
   
   void setup_mesh_worker (IntegratorExplicit<dim>&);
   
   void set_initial_condition ();
   void set_initial_condition_Qk ();
   void set_initial_condition_Pk ();
   
   std::pair<unsigned int, double> solve (LA::Vector<double> &solution, double current_residual);
   
   void compute_refinement_indicators (Vector<double> &indicator) const;
   void refine_grid (const Vector<double> &indicator);
   void refine_forward_step ();
   
   void output_results ();
   
   typedef dealii::MeshWorker::DoFInfo<dim> DoFInfo;
   typedef dealii::MeshWorker::IntegrationInfo<dim> CellInfo;
   
   // Functions for explicit integration
   void integrate_cell_term_explicit (DoFInfo& dinfo, CellInfo& info);
   void integrate_boundary_term_explicit (DoFInfo& dinfo, CellInfo& info);
   void integrate_face_term_explicit (DoFInfo& dinfo1, DoFInfo& dinfo2,
                             CellInfo& info1, CellInfo& info2);
   void assemble_system (IntegratorExplicit<dim>& integrator);

   // Functions for implicit integration
   void integrate_cell_term (DoFInfo& dinfo, CellInfo& info);
   void integrate_boundary_term (DoFInfo& dinfo, CellInfo& info);
   void integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                             CellInfo& info1, CellInfo& info2);
   void iterate_explicit (IntegratorExplicit<dim>& integrator,
                          LA::Vector<double>& newton_update,
                          double& res_norm0, double& res_norm);

   void compute_time_step ();
   void compute_time_step_cartesian ();
   void compute_time_step_q ();
   void compute_angular_momentum ();
   void compute_cell_average ();
   void apply_limiter ();
   void apply_limiter_TVB_Qk ();
   void apply_limiter_TVB_Pk ();
   void apply_limiter_minmax_Qk ();
   void apply_positivity_limiter ();
   void apply_positivity_limiter_cell
      (typename DoFHandler<dim>::active_cell_iterator& cell,
       PosLimData<dim>& data);
   void compute_shock_indicator ();
   void compute_shock_indicator_kxrcf ();
   
   void compute_mu_shock ();
   void shock_cell_term (DoFInfo& dinfo, CellInfo& info);
   void shock_boundary_term (DoFInfo& dinfo, CellInfo& info);
   void shock_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                         CellInfo& info1, CellInfo& info2);

   // The first few member variables
   // are also rather standard. Note
   // that we define a mapping
   // object to be used throughout
   // the program when assembling
   // terms (we will hand it to
   // every FEValues and
   // FEFaceValues object); the
   // mapping we use is just the
   // standard $Q_1$ mapping --
   // nothing fancy, in other words
   // -- but declaring one here and
   // using it throughout the
   // program will make it simpler
   // later on to change it if that
   // should become necessary. This
   // is, in fact, rather pertinent:
   // it is known that for
   // transsonic simulations with
   // the Euler equations,
   // computations do not converge
   // even as $h\rightarrow 0$ if
   // the boundary approximation is
   // not of sufficiently high
   // order.
   
   MPI_Comm						mpi_communicator;
   
   parallel::distributed::Triangulation<dim> triangulation;
   
   IndexSet 		      locally_owned_dofs;
   IndexSet 		      locally_relevant_dofs;
   
   const dealii::FESystem<dim>  fe;
   dealii::DoFHandler<dim>      dof_handler;
   
   // Degree zero FE for storing data on each cell
   const dealii::FE_DGQ<dim>    fe_cell;
   dealii::DoFHandler<dim>      dh_cell;

   // Iterators to neighbouring cells
   std::vector<typename dealii::DoFHandler<dim>::cell_iterator>
         lcell, rcell, bcell, tcell;
   
   // Next come a number of data
   // vectors that correspond to the
   // solution of the previous time
   // step
   // (<code>old_solution</code>),
   // the best guess of the current
   // solution
   // (<code>current_solution</code>;
   // we say <i>guess</i> because
   // the Newton iteration to
   // compute it may not have
   // converged yet, whereas
   // <code>old_solution</code>
   // refers to the fully converged
   // final result of the previous
   // time step), and a predictor
   // for the solution at the next
   // time step, computed by
   // extrapolating the current and
   // previous solution one time
   // step into the future:
   
   LA::Vector<double>			    old_solution;
   LA::Vector<double>			    current_solution;
   LA::Vector<double>			    predictor;
   LA::Vector<double>			    work1;
   LA::Vector<double>			    right_hand_side;
   LA::Vector<double>			    newton_update;
   
   std::vector< dealii::Vector<double> >	    cell_average;  
   dealii::Vector<double>       dt;
   dealii::Vector<double>       mu_shock;
   dealii::Vector<double>       shock_indicator;
   dealii::Vector<double>       jump_indicator;
   
   double                       global_dt;
   double                       elapsed_time;
   int                          time_iter;
   double                       jump_ind_min, jump_ind_max, jump_ind_avg;
   
   // This final set of member variables
   // (except for the object holding all
   // run-time parameters at the very
   // bottom and a screen output stream
   // that only prints something if
   // verbose output has been requested)
   // deals with the inteface we have in
   // this program to the Trilinos library
   // that provides us with linear
   // solvers. Similarly to including
   // PETSc matrices in step-17,
   // step-18, and step-19, all we
   // need to do is to create a Trilinos
   // sparse matrix instead of the
   // standard deal.II class. The system
   // matrix is used for the Jacobian in
   // each Newton step. Since we do not
   // intend to run this program in
   // parallel (which wouldn't be too hard
   // with Trilinos data structures,
   // though), we don't have to think
   // about anything else like
   // distributing the degrees of freedom.
   dealii::SparseMatrix<double> system_matrix;
   dealii::SparsityPattern      sparsity_pattern;

   std::vector< dealii::Vector<double> > inv_mass_matrix;
   
   Parameters::AllParameters<dim>  parameters;
   dealii::ConditionalOStream      pcout;
   TimerOutput                     computing_timer;

   // Call the appropriate numerical flux function
   template <typename InputVector>
   inline
   void numerical_normal_flux 
   (
      const dealii::Point<dim>         &normal,
      const InputVector                &Wplus,
      const InputVector                &Wminus,
      const dealii::Vector<double>     &Aplus,
      const dealii::Vector<double>     &Aminus,
      typename InputVector::value_type (&normal_flux)[EulerEquations<dim>::n_components]
   ) const
   {
      switch(parameters.flux_type)
      {
         case Parameters::Flux::lxf:
            EulerEquations<dim>::lxf_flux (normal,
                                           Wplus,
                                           Wminus,
                                           Aplus,
                                           Aminus,
                                           normal_flux);
            break;

         case Parameters::Flux::sw:
            EulerEquations<dim>::steger_warming_flux (normal,
                                                      Wplus,
                                                      Wminus,
                                                      normal_flux);
            break;

         case Parameters::Flux::kfvs:
            EulerEquations<dim>::kfvs_flux (normal,
                                            Wplus,
                                            Wminus,
                                            normal_flux);
            break;
            
         case Parameters::Flux::roe:
            EulerEquations<dim>::roe_flux (normal,
                                           Wplus,
                                           Wminus,
                                           normal_flux);
            break;
            
         case Parameters::Flux::hllc:
            EulerEquations<dim>::hllc_flux (normal,
                                            Wplus,
                                            Wminus,
                                            normal_flux);
            break;
            
         case Parameters::Flux::kep:
            EulerEquations<dim>::kep_flux (normal,
                                           Wplus,
                                           Wminus,
                                           Aplus,
                                           Aminus,
                                           normal_flux);
            break;

	      default:
            Assert (false, dealii::ExcNotImplemented());
      }
   }


   // Given a cell iterator, return the cell number
   template <typename ITERATOR>
   inline
   unsigned int cell_number (const ITERATOR &cell) const
   {
	   return cell->user_index();
   }
   
   // If cell is active, return cell average.
   // If cell is not active, return area average of child cells.
   inline
   void get_cell_average(const typename dealii::DoFHandler<dim>::cell_iterator& cell,
                         dealii::Vector<double>& avg) const
   {
      if(cell->active())
      {
         unsigned int cell_no = cell_number(cell);
         for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
            avg(c) = cell_average[cell_no][c];
      }
      else
      {  // compute average solution on child cells
         std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> child_cells =
            dealii::GridTools::get_active_child_cells< dealii::DoFHandler<dim> > (cell);
         avg = 0;
         double measure = 0;
         for(unsigned int i=0; i<child_cells.size(); ++i)
         {
            unsigned int child_cell_no = cell_number(child_cells[i]);
            for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
               avg(c) += cell_average[child_cell_no][c] * child_cells[i]->measure();
            measure += child_cells[i]->measure();
         }
         avg /= measure;
      }
   }

};

#endif
