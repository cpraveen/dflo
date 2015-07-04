#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include "claw.h"
#include "equation.h"

using namespace dealii;

// @sect4{ConservationLaw::compute_refinement_indicators}

// This function is real simple: We don't
// pretend that we know here what a good
// refinement indicator would be. Rather, we
// assume that the <code>EulerEquation</code>
// class would know about this, and so we
// simply defer to the respective function
// we've implemented there:
template <int dim>
void
ConservationLaw<dim>::
compute_refinement_indicators (dealii::Vector<double> &refinement_indicators) const 
{
   if(parameters.time_step_type == "global")
      EulerEquations<dim>::compute_refinement_indicators (dof_handler,
                                                          mapping(),
                                                          predictor,
                                                          refinement_indicators);
   else
      EulerEquations<dim>::compute_refinement_indicators (dof_handler,
                                                          mapping(),
                                                          current_solution,
                                                          refinement_indicators);
}



// @sect4{ConservationLaw::refine_grid}

// Here, we use the refinement indicators
// computed before and refine the mesh. At
// the beginning, we loop over all cells and
// mark those that we think should be
// refined:
template <int dim>
void
ConservationLaw<dim>::refine_grid (const Vector<double> &refinement_indicators)
{
   typename DoFHandler<dim>::active_cell_iterator
   cell = dof_handler.begin_active(),
   endc = dof_handler.end();
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
	   unsigned int cell_no=cell_number(cell);
	   cell->clear_coarsen_flag();
	   cell->clear_refine_flag();
	   if((cell->level() < parameters.shock_levels) &&
	      (std::fabs(refinement_indicators(cell_no)) > parameters.shock_val))
				cell->set_refine_flag();
	   else 
	      if ((cell->level() > 0) &&
			(std::fabs(refinement_indicators(cell_no)) < 0.75*parameters.shock_val))
				cell->set_coarsen_flag();
	}
   // Then we need to transfer the
   // various solution vectors from
   // the old to the new grid while we
   // do the refinement. The
   // SolutionTransfer class is our
   // friend here; it has a fairly
   // extensive documentation,
   // including examples, so we won't
   // comment much on the following
   // code. The last three lines
   // simply re-set the sizes of some
   // other vectors to the now correct
   // size:
   
   triangulation.prepare_coarsening_and_refinement();

   parallel::distributed::SolutionTransfer<dim, LA::Vector<double> > soltrans1(dof_handler), soltrans2(dof_handler);
   soltrans1.prepare_for_coarsening_and_refinement(old_solution);
   soltrans2.prepare_for_coarsening_and_refinement(predictor);
   
   triangulation.execute_coarsening_and_refinement();
   
   setup_system ();
   
   // interpolate solution to new mesh
   LA::Vector<double> distributed_solution1(locally_owned_dofs, mpi_communicator);
   LA::Vector<double> distributed_solution2(locally_owned_dofs, mpi_communicator);

   soltrans1.interpolate(distributed_solution1);
   soltrans2.interpolate(distributed_solution2);

   current_solution = distributed_solution1;
   predictor=distributed_solution2;
   
   predictor.compress(VectorOperation::insert);
   current_solution.compress(VectorOperation::insert);
}

//---------------------------------------------------------------------------
// Refine cells near the corner of forward step problem
//---------------------------------------------------------------------------
template <int dim>
void
ConservationLaw<dim>::refine_forward_step ()
{
   const double radius = 0.05;
   const Point<dim> corner(0.6, 0.2);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      cell->clear_coarsen_flag();
      cell->clear_refine_flag();
   
      Point<dim> dr = cell->center() - corner;
      if(dr.norm() < radius)
         cell->set_refine_flag();
   }
   
   triangulation.prepare_coarsening_and_refinement();
   triangulation.execute_coarsening_and_refinement ();
}

template class ConservationLaw<2>;
