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
#include "DealiiExtensions.h"

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
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;

   const FEValuesBase<dim>& fe_v    = info.fe_values();
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   const unsigned int n_q_points    = fe_v.n_quadrature_points;

   // Compute conserved variable and its gradients at all quadrature points
   Table<2,double>
      W (n_q_points, MHDEquations<dim>::n_components);
      //W (n_q_points, EulerEquations<dim>::n_components);
   
   /*Table<3,double>
      grad_W (n_q_points, EulerEquations<dim>::n_components, dim);*/

   // Compute cartesian components of flux
   //typedef double FluxMatrix[EulerEquations<dim>::n_components][dim];
   typedef double FluxMatrix[MHDEquations<dim>::n_components][dim];
   FluxMatrix *flux = new FluxMatrix[n_q_points];
   
   std::vector<Vector<double> > ext_force_values(n_q_points, Vector<double>(dim));
   parameters.external_force.vector_value_list(fe_v.get_quadrature_points(),
                                               ext_force_values);
   
   //typedef double ForcingVector[EulerEquations<dim>::n_components];
   typedef double ForcingVector[MHDEquations<dim>::n_components];
   ForcingVector *forcing = new ForcingVector[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      //for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      for(unsigned int c=0; c<MHDEquations<dim>::n_components; ++c)
      {
         W[q][c] = 0.0;
         /*for(unsigned int d=0; d<dim; ++d)
            grad_W[q][c][d] = 0.0;*/
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         
         W[q][c] += current_solution(dof_indices[i]) *
                    fe_v.shape_value_component(i, q, c);
         
         /*for (unsigned int d = 0; d < dim; ++d)
            grad_W[q][c][d] += current_solution(dof_indices[i]) *
                               fe_v.shape_grad_component(i, q, c)[d];*/
      }
      
      MHDEquations<dim>::compute_flux_matrix (W[q], flux[q]);
      MHDEquations<dim>::compute_forcing_vector (W[q], ext_force_values[q], forcing[q]);
      //EulerEquations<dim>::compute_flux_matrix (W[q], flux[q]);
      //EulerEquations<dim>::compute_forcing_vector (W[q], ext_force_values[q], forcing[q]);
   }
   
   
   // Get current cell number
   const unsigned int cell_no = cell_number (dinfo.cell);

   for (unsigned int i=0; i<dofs_per_cell; ++i)
   {
      double F_i = 0;
      
      const unsigned int
      component_i = fe_v.get_fe().system_to_component_index(i).first;
      
      for (unsigned int point=0; point<n_q_points; ++point)
      {
         for (unsigned int d=0; d<dim; d++)
            F_i -= flux[point][component_i][d] *
                   fe_v.shape_grad_component(i, point, component_i)[d] *
                   fe_v.JxW(point);
         
         // Diffusion term for shocks
         /*
         if(parameters.diffusion_coef > 0.0)
            for (unsigned int d=0; d<dim; d++)
               F_i += mu_shock(cell_no) *
                      grad_W[point][component_i][d] *
                      fe_v.shape_grad_component(i, point, component_i)[d] *
                      fe_v.JxW(point);*/
         
         F_i -= parameters.gravity *
                forcing[point][component_i] *
                fe_v.shape_value_component(i, point, component_i) *
                fe_v.JxW(point);
      }
      
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
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;
   const unsigned int& face_no = dinfo.face_number;
   const double& face_diameter = dinfo.face->diameter();
   const unsigned int& boundary_id = dinfo.face->boundary_id();
   
   const FEValuesBase<dim>& fe_v = info.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   // Get cell number of two cells adjacent to this face
   const unsigned int cell_no = cell_number (dinfo.cell);
   
   // Conservative variable value at face
   Table<2,double>
      Wplus  (n_q_points, MHDEquations<dim>::n_components),
      Wminus (n_q_points, MHDEquations<dim>::n_components);
      //Wplus  (n_q_points, EulerEquations<dim>::n_components),
      //Wminus (n_q_points, EulerEquations<dim>::n_components);

   /*Table<3,double>
      grad_Wplus  (n_q_points, EulerEquations<dim>::n_components, dim);*/
   
   // On the other side of face, we have boundary. We get Wminus from
   // boundary conditions
   
   Assert (boundary_id < Parameters::AllParameters<dim>::max_n_boundaries,
           ExcIndexRange (boundary_id, 0,
                          Parameters::AllParameters<dim>::max_n_boundaries));
   
   //typename EulerEquations<dim>::BoundaryKind boundary_kind =
   typename MHDEquations<dim>::BoundaryKind boundary_kind =
   parameters.boundary_conditions[boundary_id].kind;
   
   // Compute numerical flux at all quadrature points
   typedef double NormalFlux[MHDEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   
   /****************************************************/
   //	IN CASE OF PERIODIC BOUDNARY CONDITIONS
   /****************************************************/
   
   if(boundary_kind=MHDEquations<dim>::periodic)
   {
     // Find the neighbouring cell via map.find(key) and store the face_pair
     FacePair<dim,dim> face_pair;
     FaceCellPair<dim> cell_key(dinfo.cell, face_no);
     typename PeriodicCellMap<dim>::iterator it = periodic_map.find(cell_key);
     face_pair = it->second;
     
     // Store the info of n_cell and n_face
     typename dealii::TriaIterator<dealii::CellAccessor<dim> > n_cell;
     int face_tmp;
     if(face_pair.cell[0]==dinfo.cell)
     {
       n_cell= face_pair.cell[1];
       face_tmp = face_pair.face_idx[1];
     }
     else
     {
       n_cell = face_pair.cell[0];
       face_tmp = face_pair.face_idx[0];
     }
     const unsigned int n_face_no = face_tmp;
     const unsigned int n_cell_no = cell_number (n_cell);
     
     // Redirect from TriaAccesor to DoFCellAccessor which includes also the dof info
     typename DoFHandler<dim>::active_cell_iterator 
	n_cell_dof(&triangulation, n_cell->level(), n_cell->index(), &dof_handler);
     
     // Extracting values at the neighboring face
     QGauss<dim-1> quad_face(fe.degree+1);
     FEFaceValues<dim> fe_v_n(mapping(), fe, quad_face, update_values | update_JxW_values);
     
     std::vector<Vector<double> > n_boundary_values(n_q_points,
					   Vector<double>(MHDEquations<dim>::n_components));
     fe_v_n.reinit(n_cell_dof,n_face_no);
     fe_v_n.get_function_values (current_solution, n_boundary_values);
     const unsigned int dofs_per_cell_neighbor = fe_v_n.dofs_per_cell;

     std::vector<types::global_dof_index> dof_indices_n;
     n_cell_dof->face(n_face_no)->get_dof_indices (dof_indices_n,
						   n_cell_dof->active_fe_index());//*/
     
     // Compute fluxes at the periodic boundary
        // Compute two states on the face
   Table<2,double>
      Wplus  (n_q_points, MHDEquations<dim>::n_components),
      Wminus (n_q_points, MHDEquations<dim>::n_components);

   // Wminus is Neighbouring cell value
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      //for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      for(unsigned int c=0; c<MHDEquations<dim>::n_components; ++c)
      {
         Wplus[q][c] = 0.0;
         Wminus[q][c] = 0.0;
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][c] += current_solution(dof_indices[i]) * 
                        fe_v.shape_value_component(i, q, c);
      }
      // Check face_flip of the boundary face and change order of quadrature points if necesary
      unsigned int q1=q;
      if(face_pair.orientation[1])
	q1 = n_q_points - q;
      for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      {
	 const unsigned int c = fe_v_n.get_fe().system_to_component_index(i).first;
	 Wminus[q][c] += current_solution(dof_indices_n[i]) *
                         fe_v_n.shape_value_component(i, q1, c);
      }
      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q],
                            Wminus[q],
                            cell_average[cell_no],
                            cell_average[n_cell_no],
                            normal_fluxes[q]);//*/
   }
   
   }
   /*************************************************************/
   else
   {//*/
   
   std::vector<Vector<double> >
   //boundary_values(n_q_points, Vector<double>(EulerEquations<dim>::n_components));
   boundary_values(n_q_points, Vector<double>(MHDEquations<dim>::n_components));
   parameters.boundary_conditions[boundary_id]
   .values.vector_value_list(fe_v.get_quadrature_points(),
                             boundary_values);

   for (unsigned int q=0; q<n_q_points; ++q)
   {
      //for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      for(unsigned int c=0; c<MHDEquations<dim>::n_components; ++c)
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
      
      //EulerEquations<dim>::compute_Wminus (boundary_kind,
      MHDEquations<dim>::compute_Wminus (boundary_kind,
                                           fe_v.normal_vector(q),
                                           Wplus[q],
                                           boundary_values[q],
                                           Wminus[q]);

      // Apply bc on cell average also. This part is bit ugly.
      Table<2,double> avg (2, MHDEquations<dim>::n_components);
      //Table<2,double> avg (2, EulerEquations<dim>::n_components);
      //for(unsigned int ic=0; ic<EulerEquations<dim>::n_components; ++ic)
      for(unsigned int ic=0; ic<MHDEquations<dim>::n_components; ++ic)
      {
         avg[0][ic] = cell_average[cell_no][ic];
         avg[1][ic] = cell_average[cell_no][ic];
      }
      //EulerEquations<dim>::compute_Wminus (boundary_kind,
      MHDEquations<dim>::compute_Wminus (boundary_kind,
                                           fe_v.normal_vector(q),
                                           avg[0],
                                           boundary_values[q],
                                           avg[1]);
      Vector<double> Aplus (MHDEquations<dim>::n_components);
      Vector<double> Aminus(MHDEquations<dim>::n_components);
      //Vector<double> Aplus (EulerEquations<dim>::n_components);
      //Vector<double> Aminus(EulerEquations<dim>::n_components);
      //for(unsigned int ic=0; ic<EulerEquations<dim>::n_components; ++ic)
      for(unsigned int ic=0; ic<MHDEquations<dim>::n_components; ++ic)
      {
         Aplus[ic]  = avg[0][ic];
         Aminus[ic] = avg[1][ic];
      }

      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q],
                            Wminus[q],
                            Aplus,
                            Aminus,
                            normal_fluxes[q]);
   }
}
   
   
   // Now assemble the face term
   for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      //if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
         double F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += normal_fluxes[point][component_i] *
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

   // Get cell number of two cells adjacent to this face
   const unsigned int cell_no          = cell_number (dinfo1.cell);
   const unsigned int neighbor_cell_no = cell_number (dinfo2.cell);

   //double mu_shock_avg = 0.5 * 
      //(mu_shock(cell_no)/dinfo1.cell->diameter()
       //+
       //mu_shock(neighbor_cell_no)/dinfo2.cell->diameter());

   // Compute two states on the face
   Table<2,double>
      Wplus  (n_q_points, MHDEquations<dim>::n_components),
      Wminus (n_q_points, MHDEquations<dim>::n_components);
      //Wplus  (n_q_points, EulerEquations<dim>::n_components),
      //Wminus (n_q_points, EulerEquations<dim>::n_components);

   /*Table<3,double>
      grad_Wplus  (n_q_points, EulerEquations<dim>::n_components, dim),
      grad_Wminus (n_q_points, EulerEquations<dim>::n_components, dim);*/

   // Compute numerical flux at all quadrature points
   //typedef double NormalFlux[EulerEquations<dim>::n_components];
   typedef double NormalFlux[MHDEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   
   // Wminus is Neighbouring cell value
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      //for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      for(unsigned int c=0; c<MHDEquations<dim>::n_components; ++c)
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
   }
   
   // Now assemble the face term
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      //if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
         double F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += normal_fluxes[point][component_i] *
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
      //if (fe_v_neighbor.get_fe().has_support_on_face(i, face_no_neighbor) == true)
      {
         double F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            
            F_i -= normal_fluxes[point][component_i] *
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

}

//------------------------------------------------------------------------------
// Assemble matrices
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::assemble_system (IntegratorExplicit<dim>& integrator)
{
   TimerOutput::Scope t(computing_timer, "Assemble");

   // Compute artificial viscosity for shock capturing
   //compute_mu_shock ();

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
    
    right_hand_side.compress (VectorOperation::add);

}

template class ConservationLaw<2>;
