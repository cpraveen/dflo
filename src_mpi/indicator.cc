#include <base/quadrature_lib.h>

#include <fe/fe_values.h>

#include <dofs/dof_handler.h>

#include "equation.h"
#include "claw.h"

using namespace dealii;

//-----------------------------------------------------------------------------
// compute KXRCF shock indicator based on density or energy
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::compute_shock_indicator ()
{
   // If indicator type is "limiter" then mark all cells.
   if(parameters.shock_indicator_type == Parameters::Limiter::limiter)
   {
      shock_indicator = 1e20;
   }
   else
   {
      compute_shock_indicator_kxrcf();
   }
}

//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::compute_shock_indicator_kxrcf ()
{
   const unsigned int density_component = EulerEquations<dim>::density_component;
   const unsigned int energy_component = EulerEquations<dim>::energy_component;
   
   QGauss<dim-1> quadrature(fe.degree + 1);
   FEFaceValues<dim> fe_face_values (mapping(), fe, quadrature,
                                     update_values | update_normal_vectors);
   FEFaceValues<dim> fe_face_values_nbr (mapping(), fe, quadrature,
                                         update_values);
   FESubfaceValues<dim> fe_subface_values (mapping(), fe, quadrature,
                                           update_values | update_normal_vectors);
   FESubfaceValues<dim> fe_subface_values_nbr (mapping(), fe, quadrature,
                                               update_values);
   
   unsigned int n_q_points = quadrature.size();
   //dealii::Vector<double> face_values(n_q_points), face_values_nbr(n_q_points);
   std::vector<double> face_values(n_q_points), face_values_nbr(n_q_points);

   
   // select indicator variable
   unsigned int component;
   switch(parameters.shock_indicator_type)
   {
      case Parameters::Limiter::density:
         component = density_component;
         break;
      case Parameters::Limiter::energy:
         component = energy_component;
         break;
      default:
         component = 0;
         AssertThrow(false, ExcNotImplemented());
   }
   
   const FEValuesExtractors::Scalar variable (component);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   jump_ind_min = 1.0e20;
   jump_ind_max = 0.0;
   jump_ind_avg = 0.0;
   
   for(; cell != endc; ++cell)
   if(cell->is_locally_owned())
   {
      unsigned int c = cell_number(cell);
      double& cell_shock_ind = shock_indicator (c);
      double& cell_jump_ind = jump_indicator (c);
      
      cell_shock_ind = 0;
      cell_jump_ind = 0;
      double inflow_measure = 0;
      
      // velocity based on cell average. we use this to determine inflow/outflow
      // parts of the cell boundary.
      Point<dim> vel;
      for(unsigned int i=0; i<dim; ++i)
         vel(i) = cell_average[c][i] / cell_average[c][density_component];
      
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
         if (cell->at_boundary(f) == false)
         {
            if ((cell->neighbor(f)->level() == cell->level()) &&
                (cell->neighbor(f)->has_children() == false))
            {
               fe_face_values.reinit(cell, f);
               fe_face_values_nbr.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));
               fe_face_values[variable].get_function_values(current_solution, face_values);
               fe_face_values_nbr[variable].get_function_values(current_solution, face_values_nbr);
               for(unsigned int q=0; q<n_q_points; ++q)
               {
                  int inflow_status = (vel * fe_face_values.normal_vector(q) < 0);
                  cell_shock_ind += inflow_status *
                                    (face_values[q] - face_values_nbr[q]) *
                                    fe_face_values.JxW(q);
                  cell_jump_ind += std::pow(face_values[q] - face_values_nbr[q], 2) *
                                   fe_face_values.JxW(q);
                  inflow_measure += inflow_status * fe_face_values.JxW(q);
               }
               
            }
            else if ((cell->neighbor(f)->level() == cell->level()) &&
                     (cell->neighbor(f)->has_children() == true))
            {
               for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
               {
                  fe_subface_values.reinit (cell, f, subface);
                  fe_face_values_nbr.reinit (cell->neighbor_child_on_subface (f, subface),
                                             cell->neighbor_of_neighbor(f));
                  fe_subface_values[variable].get_function_values(current_solution, face_values);
                  fe_face_values_nbr[variable].get_function_values(current_solution, face_values_nbr);
                  for(unsigned int q=0; q<n_q_points; ++q)
                  {
                     int inflow_status = (vel * fe_subface_values.normal_vector(q) < 0);
                     cell_shock_ind += inflow_status *
                                       (face_values[q] - face_values_nbr[q]) *
                                       fe_subface_values.JxW(q);
                     cell_jump_ind += std::pow(face_values[q] - face_values_nbr[q], 2) *
                                      fe_face_values.JxW(q);
                     inflow_measure += inflow_status * fe_subface_values.JxW(q);
                  }
               }
            }
            else if (cell->neighbor_is_coarser(f))
            {
               fe_face_values.reinit(cell, f);
               fe_subface_values_nbr.reinit (cell->neighbor(f),
                                             cell->neighbor_of_coarser_neighbor(f).first,
                                             cell->neighbor_of_coarser_neighbor(f).second);
               fe_face_values[variable].get_function_values(current_solution, face_values);
               fe_subface_values_nbr[variable].get_function_values(current_solution, face_values_nbr);
               for(unsigned int q=0; q<n_q_points; ++q)
               {
                  int inflow_status = (vel * fe_face_values.normal_vector(q) < 0);
                  cell_shock_ind += inflow_status *
                                    (face_values[q] - face_values_nbr[q]) *
                                    fe_face_values.JxW(q);
                  cell_jump_ind += std::pow(face_values[q] - face_values_nbr[q], 2) *
                                   fe_face_values.JxW(q);
                  inflow_measure += inflow_status * fe_face_values.JxW(q);
               }
            }
         }
         else
         {
            // Boundary face
            // We dont do anything here since we assume solution is constant near
            // boundary.
         }
      
      // normalized shock indicator
      double cell_norm = cell_average[c][component];
      double denominator = std::pow(cell->diameter(), 0.5*(fe.degree+1)) *
                           inflow_measure *
                           cell_norm;
      cell_shock_ind = std::fabs(cell_shock_ind) / denominator;
      
      
      double dx = cell->diameter() / std::sqrt(1.0*dim);
      // TODO: normalization needs to be done properly
      cell_jump_ind = std::sqrt( cell_jump_ind / (4.0*dx) ) * cell->diameter();
      jump_ind_min = std::min(jump_ind_min, cell_jump_ind);
      jump_ind_max = std::max(jump_ind_max, cell_jump_ind);
      jump_ind_avg += cell_jump_ind;
   }
   
   jump_ind_avg /= triangulation.n_active_cells();
}

template class ConservationLaw<2>;
