#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

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
void ConservationLaw<dim>::shock_cell_term (DoFInfo& dinfo, 
                                            CellInfo& info)
{
   std::vector<unsigned int>& dof_indices = dinfo.indices;

   const FEValuesBase<dim>& fe_v    = info.fe_values();
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   const unsigned int n_q_points    = fe_v.n_quadrature_points;

   Table<2,double>
      W_theta (n_q_points, EulerEquations<dim>::n_components);

   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for (unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         W_theta[q][c] = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         
         double dof_theta = parameters.theta *
                            current_solution(dof_indices[i])
                            +
                            (1-parameters.theta) *
                            old_solution(dof_indices[i]);

         W_theta[q][c] += dof_theta *
                          fe_v.shape_value_component(i, q, c);
         
      }
   }
   
   
   // Get current cell number
   const unsigned int cell_no = cell_number (dinfo.cell);

   double density_norm = 0;
   double average_speed = 0;
   for (unsigned int point=0; point<n_q_points; ++point)
   {
      density_norm += W_theta[point][EulerEquations<dim>::density_component] *
                      fe_v.JxW(point);
      double max_speed = EulerEquations<dim>::max_eigenvalue (W_theta[point]);
      average_speed += max_speed * fe_v.JxW(point);
   }
   average_speed /= dinfo.cell->measure();

   mu_shock(cell_no) *= parameters.diffusion_coef *
                        std::pow(dinfo.cell->diameter(), 2) *
                        average_speed /
                        density_norm;

   //mu_shock(cell_no) *= parameters.diffusion_coef;

   
   //mu_shock(cell_no) = parameters.diffusion_coef *
   //                     dinfo.cell->diameter();

   /*
   mu_shock(cell_no) = parameters.diffusion_coef * average_speed * 
                       dinfo.cell->diameter();
                       */

}


//------------------------------------------------------------------------------
// Contribution from boundary faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::shock_boundary_term (DoFInfo& dinfo, 
                                                CellInfo& info)
{
   std::vector<unsigned int>& dof_indices = dinfo.indices;
   const unsigned int& face_no = dinfo.face_number;
   const unsigned int& boundary_id = dinfo.face->boundary_indicator();
   
   const FEValuesBase<dim>& fe_v = info.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   Table<2,double>
   Wplus (n_q_points, EulerEquations<dim>::n_components),
   Wminus (n_q_points, EulerEquations<dim>::n_components);

   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         Wplus[q][c] = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][component_i] += (parameters.theta *
                                   current_solution(dof_indices[i])
                                   +
                                   (1.0-parameters.theta) *
                                   old_solution(dof_indices[i])) *
                                   fe_v.shape_value_component(i, q, component_i);
      }
   }
   

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

   for (unsigned int q = 0; q < n_q_points; q++)
      EulerEquations<dim>::compute_Wminus (boundary_kind,
                                           fe_v.normal_vector(q),
                                           Wplus[q],
                                           boundary_values[q],
                                           Wminus[q]);
   
   // Compute entropy variables at quadrature points
   // We declare it as type Sacado::Fad::DFad<double> even though we dont need its
   // derivative. Otherwise there is error in assignment due to different types.
   typedef double EntropyVar[EulerEquations<dim>::n_components];
   EntropyVar *Vplus = new EntropyVar[n_q_points];
   EntropyVar *Vminus= new EntropyVar[n_q_points];
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      EulerEquations<dim>::entropy_var (Wplus [q], Vplus [q]);
      EulerEquations<dim>::entropy_var (Wminus[q], Vminus[q]);
   }

   // Compute integral of shock indicator on face
   double jump = 0;
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      double ds = 0;
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         ds += (Vplus[q][c] - Vminus[q][c]) * (Wplus[q][c] - Wminus[q][c]);
      jump += ds * fe_v.JxW(q);
   }

   // Add shock indicator to the two cells
   mu_shock(cell_number(dinfo.cell)) += jump;

   delete[] Vplus;
   delete[] Vminus;
}



//------------------------------------------------------------------------------
// Contribution from interior faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::shock_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                                            CellInfo& info1, CellInfo& info2)
{
   std::vector<unsigned int>& dof_indices = dinfo1.indices;
   const unsigned int& face_no = dinfo1.face_number;
   
   std::vector<unsigned int>& dof_indices_neighbor = dinfo2.indices;
   const unsigned int& face_no_neighbor = dinfo2.face_number;

   const FEValuesBase<dim>& fe_v = info1.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   const FEValuesBase<dim>& fe_v_neighbor = info2.fe_values();
   const unsigned int dofs_per_cell_neighbor = fe_v_neighbor.dofs_per_cell;

   Table<2,double>
   Wplus (n_q_points, EulerEquations<dim>::n_components),
   Wminus (n_q_points, EulerEquations<dim>::n_components);
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         Wplus[q][c] = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][component_i] += (parameters.theta *
                                   current_solution(dof_indices[i])
                                   +
                                   (1.0-parameters.theta) *
                                   old_solution(dof_indices[i])) *
                                   fe_v.shape_value_component(i, q, component_i);
      }
   }
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         Wminus[q][c] = 0;
      for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      {
         const unsigned int component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
         Wminus[q][component_i] += (parameters.theta *
                                    current_solution(dof_indices_neighbor[i])
                                    +
                                    (1.0-parameters.theta) *
                                    old_solution(dof_indices_neighbor[i]))*
                                    fe_v_neighbor.shape_value_component(i, q, component_i);
      }
   }
   
   
   // Compute entropy variables at quadrature points
   // We declare it as type Sacado::Fad::DFad<double> even though we dont need its
   // derivative. Otherwise there is error in assignment due to different types.
   typedef double EntropyVar[EulerEquations<dim>::n_components];
   EntropyVar *Vplus = new EntropyVar[n_q_points];
   EntropyVar *Vminus= new EntropyVar[n_q_points];
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      EulerEquations<dim>::entropy_var (Wplus [q], Vplus [q]);
      EulerEquations<dim>::entropy_var (Wminus[q], Vminus[q]);
   }

   // Compute integral of shock indicator on face
   double jump = 0;
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      double ds = 0;
      for(unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
         ds += (Vplus[q][c] - Vminus[q][c]) * (Wplus[q][c] - Wminus[q][c]);
      jump += ds * fe_v.JxW(q);
   }

   // Add shock indicator to the two cells
   mu_shock(cell_number(dinfo1.cell)) += jump;
   mu_shock(cell_number(dinfo2.cell)) += jump;

   delete[] Vplus;
   delete[] Vminus;

}

//------------------------------------------------------------------------------
// Assemble matrices
//------------------------------------------------------------------------------
/*template <int dim>
void ConservationLaw<dim>::compute_mu_shock ()
{
   if(parameters.diffusion_coef == 0.0)
   {
      mu_shock = 0;
      return;
   }

   IntegratorImplicit<dim> integrator(dof_handler);
   const unsigned int n_gauss_points = fe.degree + 1;
   integrator.info_box.initialize_gauss_quadrature(n_gauss_points,
                                                   n_gauss_points,
                                                   n_gauss_points);
   
   integrator.info_box.initialize_update_flags ();
   integrator.info_box.add_update_flags_all (update_values | 
                                             update_quadrature_points |
                                             update_JxW_values);
   //integrator.info_box.add_update_flags_cell     (update_gradients);
   //integrator.info_box.add_update_flags_boundary (update_normal_vectors);
   //integrator.info_box.add_update_flags_face     (update_normal_vectors);
   
   integrator.info_box.initialize (fe, mapping());
   
   integrator.assembler.initialize (system_matrix, right_hand_side);

   mu_shock = 0;

   MeshWorker::loop<dim,dim,MeshWorker::DoFInfo<dim>,MeshWorker::IntegrationInfoBox<dim> >
   (dof_handler.begin_active(),
    dof_handler.end(),
    integrator.dof_info, 
    integrator.info_box,
    0,
    boost::bind(&ConservationLaw<dim>::shock_boundary_term,
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::shock_face_term,
                this, _1, _2, _3, _4),
    integrator.assembler);

   MeshWorker::loop<dim,dim,MeshWorker::DoFInfo<dim>,MeshWorker::IntegrationInfoBox<dim> >
   (dof_handler.begin_active(),
    dof_handler.end(),
    integrator.dof_info, 
    integrator.info_box,
    boost::bind(&ConservationLaw<dim>::shock_cell_term, 
                this, _1, _2),
    0,
    0,
    integrator.assembler);

   return;

   // Smooth viscosity by averaging
   FE_Q<dim> fe1 (1);
   DoFHandler<dim> dh1(triangulation);
   dh1.distribute_dofs (fe1);
   Vector<double> mu_tmp (dh1.n_dofs());
   mu_tmp = 0;
   std::vector<unsigned int> dof_indices (fe1.dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
     cell = dh1.begin_active(),
     endc = dh1.end();
   for (; cell!=endc; ++cell)
     {
         const unsigned int cell_no = cell_number (cell);

         cell->get_dof_indices (dof_indices);
         for(unsigned int i=0; i<fe1.dofs_per_cell; ++i)
            mu_tmp(dof_indices[i]) += 0.25 * mu_shock (cell_no);
     }

   cell = dh1.begin_active();
   endc = dh1.end();
   for (; cell!=endc; ++cell)
     {
         const unsigned int cell_no = cell_number (cell);
         mu_shock (cell_no) = 0;

         cell->get_dof_indices (dof_indices);
         for(unsigned int i=0; i<fe1.dofs_per_cell; ++i)
            mu_shock(cell_no) += 0.25 * mu_tmp(dof_indices[i]);
     }
}//*/

template class ConservationLaw<2>;
