#ifndef __CLAW_H__
#define __CLAW_H__

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>

#include <lac/vector.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/precondition_block.h>

#include <grid/tria.h>
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

#include <vector>

#include "parameters.h"
#include "integrator.h"

template <int dim>
struct CellData
{
   typename dealii::DoFHandler<dim>::cell_iterator lcell, rcell, bcell, tcell;
   unsigned int number;
   typename EulerEquations<dim>::BoundaryKind lbc, rbc, bbc, tbc;
};

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
   const dealii::Mapping<dim,dim>& mapping() const;
   void compute_cartesian_mesh_size ();
   void compute_inv_mass_matrix();
   void setup_system ();
   
   void setup_mesh_worker (IntegratorImplicit<dim>&);
   void setup_mesh_worker (IntegratorExplicit<dim>&);
   
   void set_initial_condition ();
   void set_initial_condition_Qk ();
   void set_initial_condition_Pk ();
   
   std::pair<unsigned int, double> solve (dealii::Vector<double> &solution, double current_residual);
   
   void compute_refinement_indicators (dealii::Vector<double> &indicator) const;
   void refine_grid (const dealii::Vector<double> &indicator);
   void refine_forward_step ();
   
   void output_results () const;
   
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
   void assemble_system (IntegratorImplicit<dim>& integrator);
   void iterate_explicit (IntegratorExplicit<dim>& integrator,
                          dealii::Vector<double>& newton_update,
                          double& res_norm0, double& res_norm);
   void iterate_mood     (IntegratorExplicit<dim>& integrator,
                          dealii::Vector<double>& newton_update,
                          double& res_norm0, double& res_norm);
   void iterate_implicit (IntegratorImplicit<dim>& integrator,
                          dealii::Vector<double>& newton_update,
                          double& res_norm0, double& res_norm);

   void compute_time_step ();
   void compute_time_step_cartesian ();
   void compute_time_step_q ();
   void compute_angular_momentum ();
   void compute_cell_average ();
   void apply_limiter ();
   void apply_limiter_TVB_Qk ();
   void apply_limiter_minmax_Qk ();
   void apply_limiter_TVB_Pk ();
   void apply_positivity_limiter ();
   void compute_shock_indicator ();
   void compute_shock_indicator_u2 ();
   void compute_shock_indicator_kxrcf ();
   
   void compute_reduction_matrices();
   void compute_min_max_mood_var();
   bool apply_mood(unsigned int&, unsigned int&, unsigned int&);
   void reduce_degree(const typename dealii::DoFHandler<dim>::cell_iterator&,
                      const unsigned int,
                      dealii::FEValues<dim>&);
   void reduce_degree_Pk(const typename dealii::DoFHandler<dim>::cell_iterator&,
                         const unsigned int,
                         dealii::FEValues<dim>&);
   void reduce_degree_Qk(const typename dealii::DoFHandler<dim>::cell_iterator&,
                         const unsigned int,
                         dealii::FEValues<dim>&);
   void get_mood_second_derivatives(const typename dealii::DoFHandler<dim>::cell_iterator &cell,
                                    std::vector<double>& D2);
   bool test_u2(const typename dealii::DoFHandler<dim>::cell_iterator &cell);
   
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
   dealii::Triangulation<dim>   triangulation;
   
   const dealii::FESystem<dim>  fe;
   dealii::DoFHandler<dim>      dof_handler;
   
   // Degree zero FE for storing data on each cell
   const dealii::FE_DGQ<dim>    fe_cell;
   dealii::DoFHandler<dim>      dh_cell;

   // Iterators to neighbouring cells
   std::vector<CellData<dim> > cell_data;
   
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
   dealii::Vector<double>       old_solution;
   dealii::Vector<double>       current_solution;
   dealii::Vector<double>       predictor;
   dealii::Vector<double>       work1;
   std::vector< dealii::Vector<double> >       cell_average;
   
   dealii::Vector<double>       right_hand_side;

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
   dealii::PreconditionBlockJacobi<dealii::SparseMatrix<double>,double> inv_mass_matrix_full;
   
   // MOOD data
   std::vector<unsigned int> cell_degree;
   std::vector<bool> re_update;
   std::vector<dealii::FullMatrix<double> > Rmatrix;
   dealii::Vector<double> min_mood_var, max_mood_var;

   std::vector< dealii::Vector<double> > inv_mass_matrix_diag;

   // For FE_DGP, maps dof index to its degree
   std::vector<unsigned int> index_to_degree;
   
   Parameters::AllParameters<dim>  parameters;
   dealii::ConditionalOStream      verbose_cout;

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
         auto child_cells =
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
