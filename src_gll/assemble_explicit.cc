#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_dgq.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "claw.h"

using namespace dealii;

//------------------------------------------------------------------------------
// Contribution of volume integral terms
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_cell_term_explicit
(
   DoFInfo  &dinfo, 
   CellInfo &info
)
{
   // Get current cell number
   const unsigned int cell_no = cell_number (dinfo.cell);
   if(re_update[cell_no] == false) return;
   
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;

   const FEValuesBase<dim>& fe_v    = info.fe_values();
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   const unsigned int n_q_points    = fe_v.n_quadrature_points;

   // Compute conserved variable and its gradients at all quadrature points
   Table<2,double>
      W (n_q_points, EulerEquations<dim>::n_components);

   // Compute cartesian components of flux
   typedef double FluxMatrix[EulerEquations<dim>::n_components][dim];
   FluxMatrix *flux = new FluxMatrix[n_q_points];
   
   typedef double ForcingVector[EulerEquations<dim>::n_components];
   ForcingVector *forcing = new ForcingVector[n_q_points];
   std::vector<double> potential(n_q_points), exp_potential(n_q_points);
   std::vector< Vector<double> > grad_exp(n_q_points, Vector<double>(dim));
   
   parameters.potential.value_list(fe_v.get_quadrature_points(),
                                   potential);
   double rho_avg = cell_average[cell_no][EulerEquations<dim>::density_component];
   double pre_avg = EulerEquations<dim>::template compute_pressure<double>(cell_average[cell_no]);
   double RT = pre_avg / rho_avg;
   
   for (unsigned int i=0; i<dofs_per_cell; ++i)
   {
      const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
      const unsigned int q = fe_v.get_fe().system_to_component_index(i).second;
      
      W[q][c] = current_solution(dof_indices[i]);
   }
      
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      EulerEquations<dim>::compute_flux_matrix (W[q], flux[q]);
      //EulerEquations<dim>::compute_forcing_vector (W[q], forcing[q]);
      exp_potential[q] = std::exp(-potential[q]/RT);
   }
   
   // compute gradient of exp(-phi/RT)
   // we use shape_grad corresponding to component = 0
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int d=0; d<dim; ++d) grad_exp[q][d] = 0;
      for(unsigned int i=0; i<n_q_points; ++i)
      {
         const unsigned int index = fe_v.get_fe().component_to_system_index(0, i);
         for(unsigned int d=0; d<dim; ++d)
            grad_exp[q][d] += exp_potential[i] *
                              fe_v.shape_grad_component(index, q, 0)[d];
      }
      
      double pressure = W[q][EulerEquations<dim>::density_component] * RT;
      double factor = pressure / exp_potential[q];
      
      forcing[q][EulerEquations<dim>::energy_component] = 0.0;
      for(unsigned int d=0; d<dim; ++d)
      {
         forcing[q][d] = factor * grad_exp[q][d];
         double u = W[q][d] / W[q][EulerEquations<dim>::density_component];
         forcing[q][EulerEquations<dim>::energy_component] += u * forcing[q][d];
      }
      forcing[q][EulerEquations<dim>::density_component] = 0.0;
   }
   
   for (unsigned int i=0; i<dofs_per_cell; ++i)
   {
      double F_i = 0;
      
      const unsigned int
      component_i = fe_v.get_fe().system_to_component_index(i).first;
      const unsigned int
      support_i = fe_v.get_fe().system_to_component_index(i).second;

      for (unsigned int point=0; point<n_q_points; ++point)
      {
         const unsigned int index = fe_v.get_fe().component_to_system_index(component_i, point);
         for (unsigned int d=0; d<dim; d++)
            F_i += flux[point][component_i][d] *
                   fe_v.shape_grad_component(index, support_i, component_i)[d];
      }

      F_i = (F_i - forcing[support_i][component_i]) * fe_v.JxW(support_i);
      
      local_vector (i) -= F_i;
   }
   
   delete[] forcing;
   delete[] flux;
   
}

//------------------------------------------------------------------------------
// Contribution of volume integral terms
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_cell_term_explicit_polytropic
(
 DoFInfo  &dinfo,
 CellInfo &info
 )
{
   // Polytropic index
   const double nu = 1.2;
   
   // Get current cell number
   const unsigned int cell_no = cell_number (dinfo.cell);
   if(re_update[cell_no] == false) return;
   
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;
   
   const FEValuesBase<dim>& fe_v    = info.fe_values();
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   const unsigned int n_q_points    = fe_v.n_quadrature_points;
   
   // Compute conserved variable and its gradients at all quadrature points
   Table<2,double>
   W (n_q_points, EulerEquations<dim>::n_components);
   
   // Compute cartesian components of flux
   typedef double FluxMatrix[EulerEquations<dim>::n_components][dim];
   FluxMatrix *flux = new FluxMatrix[n_q_points];
   
   typedef double ForcingVector[EulerEquations<dim>::n_components];
   ForcingVector *forcing = new ForcingVector[n_q_points];
   std::vector<double> potential(n_q_points), exp_potential(n_q_points);
   std::vector< Vector<double> > grad_exp(n_q_points, Vector<double>(dim));
   
   parameters.potential.value_list(fe_v.get_quadrature_points(),
                                   potential);
   
   for (unsigned int i=0; i<dofs_per_cell; ++i)
   {
      const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
      const unsigned int q = fe_v.get_fe().system_to_component_index(i).second;
      
      W[q][c] = current_solution(dof_indices[i]);
   }
   
   double alpha, beta = -1.0e20;
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      EulerEquations<dim>::compute_flux_matrix (W[q], flux[q]);
      //EulerEquations<dim>::compute_forcing_vector (W[q], forcing[q]);
      
      double pressure = EulerEquations<dim>::template compute_pressure<double>(W[q]);
      double ff = nu/(nu - 1.0) * pressure/W[q][EulerEquations<dim>::density_component] + potential[q];
      if(ff > beta)
      {
         beta = ff;
         alpha = pressure / std::pow(W[q][EulerEquations<dim>::density_component], nu);
      }
   }
   
   for (unsigned int q=0; q<n_q_points; ++q)
      exp_potential[q] = std::pow((nu-1.0)*(beta - potential[q])/(nu*alpha), nu/(nu-1.0));
   
   // compute gradient of exp(-phi/RT)
   // we use shape_grad corresponding to component = 0
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int d=0; d<dim; ++d) grad_exp[q][d] = 0;
      for(unsigned int i=0; i<n_q_points; ++i)
      {
         const unsigned int index = fe_v.get_fe().component_to_system_index(0, i);
         for(unsigned int d=0; d<dim; ++d)
            grad_exp[q][d] += exp_potential[i] *
                              fe_v.shape_grad_component(index, q, 0)[d];
      }
      
      double factor = (nu-1.0)/nu * W[q][EulerEquations<dim>::density_component]
                    * (beta - potential[q]) / exp_potential[q];
      
      forcing[q][EulerEquations<dim>::energy_component] = 0.0;
      for(unsigned int d=0; d<dim; ++d)
      {
         forcing[q][d] = factor * grad_exp[q][d];
         double u = W[q][d] / W[q][EulerEquations<dim>::density_component];
         forcing[q][EulerEquations<dim>::energy_component] += u * forcing[q][d];
      }
      forcing[q][EulerEquations<dim>::density_component] = 0.0;
   }
   
   for (unsigned int i=0; i<dofs_per_cell; ++i)
   {
      double F_i = 0;
      
      const unsigned int
      component_i = fe_v.get_fe().system_to_component_index(i).first;
      const unsigned int
      support_i = fe_v.get_fe().system_to_component_index(i).second;
      
      for (unsigned int point=0; point<n_q_points; ++point)
      {
         const unsigned int index = fe_v.get_fe().component_to_system_index(component_i, point);
         for (unsigned int d=0; d<dim; d++)
            F_i += flux[point][component_i][d] *
                   fe_v.shape_grad_component(index, support_i, component_i)[d];
      }
      
      F_i = (F_i - forcing[support_i][component_i]) * fe_v.JxW(support_i);
      
      local_vector (i) -= F_i;
   }
   
   delete[] forcing;
   delete[] flux;
   
}

//------------------------------------------------------------------------------
// Contribution from boundary faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_boundary_term_explicit
(
   DoFInfo  &dinfo, 
   CellInfo &info
)
{
   // Get cell number of two cells adjacent to this face
   const unsigned int cell_no = cell_number (dinfo.cell);
   if(re_update[cell_no] == false) return;
   
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;
   const unsigned int& face_no = dinfo.face_number;
   const double& face_diameter = dinfo.face->diameter();
   const unsigned int& boundary_id = dinfo.face->boundary_id();
   
   const FEValuesBase<dim>& fe_v = info.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   // Conservative variable value at face
   Table<2,double>
      Wplus  (n_q_points, EulerEquations<dim>::n_components),
      Wminus (n_q_points, EulerEquations<dim>::n_components);

   /*Table<3,double>
      grad_Wplus  (n_q_points, EulerEquations<dim>::n_components, dim);*/
   
   // On the other side of face, we have boundary. We get Wminus from
   // boundary conditions
   
   Assert (boundary_id < Parameters::AllParameters<dim>::max_n_boundaries,
           ExcIndexRange (boundary_id, 0,
                          Parameters::AllParameters<dim>::max_n_boundaries));
   
   std::vector<Vector<double> >
   boundary_values(n_q_points, Vector<double>(EulerEquations<dim>::n_components));
   parameters.boundary_conditions[boundary_id]
   .values.vector_value_list(fe_v.get_quadrature_points(),
                             boundary_values);
   
   
   typename EulerEquations<dim>::BoundaryKind boundary_kind =
   parameters.boundary_conditions[boundary_id].kind;
   
   
   // Compute numerical flux at all quadrature points
   typedef double NormalFlux[EulerEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   NormalFlux *numerical_normal_fluxes = new NormalFlux[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      {
         Wplus[q][c] = 0.0;
         /*for(unsigned int d=0; d<dim; ++d)
            grad_Wplus[q][c][d] = 0.0;*/
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][c] += current_solution(dof_indices[i]) *
                        fe_v.shape_value_component(i, q, c);

         /*for (unsigned int d=0; d<dim; ++d)
            grad_Wplus[q][c][d] += current_solution(dof_indices[i]) *
                                   fe_v.shape_grad_component(i, q, c)[d];*/
      }
      
      EulerEquations<dim>::compute_Wminus (boundary_kind,
                                           fe_v.normal_vector(q),
                                           Wplus[q],
                                           boundary_values[q],
                                           Wminus[q]);
      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q],
                            Wminus[q],
                            cell_average[cell_no],
                            cell_average[cell_no],
                            numerical_normal_fluxes[q]);
      
      EulerEquations<dim>::normal_flux(Wplus[q],
                                       fe_v.normal_vector(q),
                                       normal_fluxes[q]);
   }
   
   // Now assemble the face term
   for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
         double F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += (numerical_normal_fluxes[point][component_i] - normal_fluxes[point][component_i]) *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);

            /*
            F_i += 5.0 * mu_shock(cell_no) / dinfo.cell->diameter() *
                   (Wplus[point][component_i] - Wminus[point][component_i]) *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);

            for(unsigned int d=0; d<dim; ++d)
               F_i -= mu_shock(cell_no)         * grad_Wplus[point][component_i][d] *
                     fe_v.normal_vector(point)[d] *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point)
                     +
                     mu_shock(cell_no) * 
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.normal_vector(point)[d] *
                     (Wplus[point][component_i] - Wminus[point][component_i]) *
                     fe_v.JxW(point);
                     */
         }
         
         local_vector (i) -= F_i;
      }
   
   delete[] normal_fluxes;
   delete[] numerical_normal_fluxes;

}



//------------------------------------------------------------------------------
// Contribution from interior faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_face_term_explicit
(
   DoFInfo  &dinfo1, 
   DoFInfo  &dinfo2,
   CellInfo &info1, 
   CellInfo &info2
)
{
   // Get cell number of two cells adjacent to this face
   const unsigned int cell_no          = cell_number (dinfo1.cell);
   const unsigned int neighbor_cell_no = cell_number (dinfo2.cell);
   if(re_update[cell_no] == false && re_update[neighbor_cell_no] == false) return;
   
   Vector<double>& local_vector   = dinfo1.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo1.indices;
   const unsigned int& face_no = dinfo1.face_number;
   const double& face_diameter = dinfo1.face->diameter();
   
   Vector<double>& local_vector_neighbor   = dinfo2.vector(0).block(0);
   std::vector<unsigned int>& dof_indices_neighbor = dinfo2.indices;
   const unsigned int& face_no_neighbor = dinfo2.face_number;

   const FEValuesBase<dim>& fe_v = info1.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   const FEValuesBase<dim>& fe_v_neighbor = info2.fe_values();
   const unsigned int dofs_per_cell_neighbor = fe_v_neighbor.dofs_per_cell;

   double mu_shock_avg = 0.5 * 
      (mu_shock(cell_no)/dinfo1.cell->diameter()
       +
       mu_shock(neighbor_cell_no)/dinfo2.cell->diameter());

   // Compute two states on the face
   Table<2,double>
      Wplus  (n_q_points, EulerEquations<dim>::n_components),
      Wminus (n_q_points, EulerEquations<dim>::n_components);

   /*Table<3,double>
      grad_Wplus  (n_q_points, EulerEquations<dim>::n_components, dim),
      grad_Wminus (n_q_points, EulerEquations<dim>::n_components, dim);*/

   // Compute numerical flux at all quadrature points
   typedef double NormalFlux[EulerEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   NormalFlux *normal_fluxes_plus = new NormalFlux[n_q_points];
   NormalFlux *normal_fluxes_minus = new NormalFlux[n_q_points];

   // Wminus is Neighbouring cell value
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      {
         Wplus[q][c] = 0.0;
         /*for(unsigned int d=0; d<dim; ++d)
            grad_Wplus[q][c][d] = 0.0;*/
         
         Wminus[q][c] = 0.0;
         /*for(unsigned int d=0; d<dim; ++d)
          grad_Wminus[q][c][d] = 0.0;*/
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][c] += current_solution(dof_indices[i]) * 
                        fe_v.shape_value_component(i, q, c);

         /*for (unsigned int d = 0; d < dim; d++)
            grad_Wplus[q][c][d] += current_solution(dof_indices[i]) *
                                   fe_v.shape_grad_component(i, q, c)[d];*/
      }
      for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      {
         const unsigned int c = fe_v_neighbor.get_fe().system_to_component_index(i).first;
         Wminus[q][c] += current_solution(dof_indices_neighbor[i]) *
                         fe_v_neighbor.shape_value_component(i, q, c);
         
         /*for (unsigned int d = 0; d < dim; d++)
          grad_Wminus[q][c][d] += current_solution(dof_indices_neighbor[i]) *
                                  fe_v_neighbor.shape_grad_component(i, q, c)[d];*/
      }
      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q],
                            Wminus[q],
                            cell_average[cell_no],
                            cell_average[neighbor_cell_no],
                            normal_fluxes[q]);
      EulerEquations<dim>::normal_flux(Wplus[q],
                                       fe_v.normal_vector(q),
                                       normal_fluxes_plus[q]);
      EulerEquations<dim>::normal_flux(Wminus[q],
                                       fe_v.normal_vector(q),
                                       normal_fluxes_minus[q]);
   }
   
   // Now assemble the face term
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
         double F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += (normal_fluxes[point][component_i] - normal_fluxes_plus[point][component_i]) *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);

            /*
            F_i += 5.0 * mu_shock_avg *
                   (Wplus[point][component_i] - Wminus[point][component_i]) *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);

            for(unsigned int d=0; d<dim; ++d)
               F_i -= 0.5 *
                     (mu_shock(cell_no)         * grad_Wplus[point][component_i][d] +
                     mu_shock(neighbor_cell_no) * grad_Wminus[point][component_i][d]) *
                     fe_v.normal_vector(point)[d] *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point)
                     +
                     0.5 *
                     mu_shock(cell_no) * 
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.normal_vector(point)[d] *
                     (Wplus[point][component_i] - Wminus[point][component_i]) *
                     fe_v.JxW(point);
                     */
         }
         
         local_vector(i) -= F_i;
      }
   
   // Contributions to neighbouring cell
   for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      if (fe_v_neighbor.get_fe().has_support_on_face(i, face_no_neighbor) == true)
      {
         double F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            
            F_i -= (normal_fluxes[point][component_i] - normal_fluxes_minus[point][component_i]) *
                   fe_v_neighbor.shape_value_component(i, point, component_i) *
                   fe_v_neighbor.JxW(point);

            /*
            F_i -= 5.0 * mu_shock_avg *
                   (Wplus[point][component_i] - Wminus[point][component_i]) *
                   fe_v_neighbor.shape_value_component(i, point, component_i) *
                   fe_v_neighbor.JxW(point);

            for(unsigned int d=0; d<dim; ++d)
               F_i += 0.5 *
                     (mu_shock(cell_no)         * grad_Wplus[point][component_i][d] +
                     mu_shock(neighbor_cell_no) * grad_Wminus[point][component_i][d]) *
                     fe_v.normal_vector(point)[d] *
                     fe_v_neighbor.shape_value_component(i, point, component_i) *
                     fe_v_neighbor.JxW(point)
                     -
                     0.5 *
                     mu_shock(neighbor_cell_no) * 
                     fe_v_neighbor.shape_grad_component(i, point, component_i)[d] *
                     fe_v.normal_vector(point)[d] *
                     (Wplus[point][component_i] - Wminus[point][component_i]) *
                     fe_v_neighbor.JxW(point);
                     */
         }
         
         local_vector_neighbor (i) -= F_i;
      }
   
   delete[] normal_fluxes;
   delete[] normal_fluxes_plus;
   delete[] normal_fluxes_minus;

}

//------------------------------------------------------------------------------
// Assemble matrices
// For polytropic scheme, pass
//        ConservationLaw<dim>::integrate_cell_term_explicit_polytropic
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::assemble_system (IntegratorExplicit<dim>& integrator)
{
   // Compute artificial viscosity for shock capturing
   compute_mu_shock ();

   right_hand_side = 0;

   MeshWorker::loop<dim,dim,MeshWorker::DoFInfo<dim>,MeshWorker::IntegrationInfoBox<dim> >
   (dof_handler.begin_active(),
    dof_handler.end(),
    integrator.dof_info, 
    integrator.info_box,
    boost::bind(&ConservationLaw<dim>::integrate_cell_term_explicit,
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::integrate_boundary_term_explicit,
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::integrate_face_term_explicit,
                this, _1, _2, _3, _4),
    integrator.assembler);
}

template class ConservationLaw<2>;
