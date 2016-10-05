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

#include <Sacado.hpp>

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
void ConservationLaw<dim>::integrate_cell_term (DoFInfo& dinfo, 
                                                CellInfo& info)
{
   FullMatrix<double>& local_matrix = dinfo.matrix(0).matrix;
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;

   const FEValuesBase<dim>& fe_v    = info.fe_values();
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   const unsigned int n_q_points    = fe_v.n_quadrature_points;

   // Get time step for current cell
   const double time_step = dt( cell_number(dinfo.cell) );
   
   Table<2,Sacado::Fad::DFad<double> >
      dWdt (n_q_points, EulerEquations<dim>::n_components);
   
   Table<2,Sacado::Fad::DFad<double> >
      W_theta (n_q_points, EulerEquations<dim>::n_components);
   
   Table<3,Sacado::Fad::DFad<double> >
      grad_W (n_q_points, EulerEquations<dim>::n_components, dim);

   // Next, we have to define the independent
   // variables that we will try to determine
   // by solving a Newton step. These
   // independent variables are the values of
   // the local degrees of freedom which we
   // extract here:
   std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
   
   // The next step incorporates all the
   // magic: we declare a subset of the
   // autodifferentiation variables as
   // independent degrees of freedom, whereas
   // all the other ones remain dependent
   // functions. These are precisely the local
   // degrees of freedom just extracted. All
   // calculations that reference them (either
   // directly or indirectly) will accumulate
   // sensitivies with respect to these
   // variables.
   //
   // In order to mark the variables as
   // independent, the following does the
   // trick, marking
   // <code>independent_local_dof_values[i]</code>
   // as the $i$th independent variable out of
   // a total of <code>dofs_per_cell</code>:
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i].diff (i, dofs_per_cell);
   
   // After all these declarations, let us
   // actually compute something. First, the
   // values of <code>W</code>,
   // <code>W_old</code>,
   // <code>W_theta</code>, and
   // <code>grad_W</code>, which we can
   // compute from the local DoF values by
   // using the formula $W(x_q)=\sum_i \mathbf
   // W_i \Phi_i(x_q)$, where $\mathbf W_i$ is
   // the $i$th entry of the (local part of
   // the) solution vector, and $\Phi_i(x_q)$
   // the value of the $i$th vector-valued
   // shape function evaluated at quadrature
   // point $x_q$. The gradient can be
   // computed in a similar way.
   //
   // Ideally, we could compute this
   // information using a call into something
   // like FEValues::get_function_values and
   // FEValues::get_function_grads, but since
   // (i) we would have to extend the FEValues
   // class for this, and (ii) we don't want
   // to make the entire
   // <code>old_solution</code> vector fad
   // types, only the local cell variables, we
   // explicitly code the loop above. Before
   // this, we add another loop that
   // initializes all the fad variables to
   // zero:
   
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         
         dWdt[q][c] += (independent_local_dof_values[i] -
                        old_solution(dof_indices[i])) / time_step *
                        fe_v.shape_value_component(i, q, c);

         Sacado::Fad::DFad<double> dof_theta = parameters.theta *
                                               independent_local_dof_values[i]
                                             +
                                               (1-parameters.theta) *
                                               old_solution(dof_indices[i]);

         W_theta[q][c] += dof_theta *
                          fe_v.shape_value_component(i, q, c);
         
         for (unsigned int d = 0; d < dim; d++)
            grad_W[q][c][d] += dof_theta *
                               fe_v.shape_grad_component(i, q, c)[d];
      }
   
   
   // Next, in order to compute the cell
   // contributions, we need to evaluate
   // $F(\tilde{\mathbf w})$ and
   // $G(\tilde{\mathbf w})$ at all quadrature
   // points. To store these, we also need to
   // allocate a bit of memory. Note that we
   // compute the flux matrices and right hand
   // sides in terms of autodifferentiation
   // variables, so that the Jacobian
   // contributions can later easily be
   // computed from it:
   typedef Sacado::Fad::DFad<double> FluxMatrix[EulerEquations<dim>::n_components][dim];
   FluxMatrix *flux = new FluxMatrix[n_q_points];
   
   typedef Sacado::Fad::DFad<double> ForcingVector[EulerEquations<dim>::n_components];
   ForcingVector *forcing = new ForcingVector[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      EulerEquations<dim>::compute_flux_matrix (W_theta[q], flux[q]);
      EulerEquations<dim>::compute_forcing_vector (W_theta[q], forcing[q]);
   }
   
   
   // Get current cell number
   const unsigned int cell_no = cell_number (dinfo.cell);

   // We now have all of the pieces in place,
   // so perform the assembly.  We have an
   // outer loop through the components of the
   // system, and an inner loop over the
   // quadrature points, where we accumulate
   // contributions to the $i$th residual
   // $F_i$. The general formula for this
   // residual is given in the introduction
   // and at the top of this function. We can,
   // however, simplify it a bit taking into
   // account that the $i$th (vector-valued)
   // test function $\mathbf{z}_i$ has in
   // reality only a single nonzero component
   // (more on this topic can be found in the
   // @ref vector_valued module). It will be
   // represented by the variable
   // <code>component_i</code> below. With
   // this, the residual term can be
   // re-written as $F_i =
   // \left(\frac{(\mathbf{w}_{n+1} -
   // \mathbf{w}_n)_{\text{component\_i}}}{\delta
   // t},(\mathbf{z}_i)_{\text{component\_i}}\right)_K$
   // $- \sum_{d=1}^{\text{dim}}
   // \left(\mathbf{F}
   // (\tilde{\mathbf{w}})_{\text{component\_i},d},
   // \frac{\partial(\mathbf{z}_i)_{\text{component\_i}}}
   // {\partial x_d}\right)_K$ $+
   // \sum_{d=1}^{\text{dim}} h^{\eta}
   // \left(\frac{\partial
   // \mathbf{w}_{\text{component\_i}}}{\partial
   // x_d} , \frac{\partial
   // (\mathbf{z}_i)_{\text{component\_i}}}{\partial
   // x_d} \right)_K$
   // $-(\mathbf{G}(\tilde{\mathbf{w}}
   // )_{\text{component\_i}},
   // (\mathbf{z}_i)_{\text{component\_i}})_K$,
   // where integrals are understood to be
   // evaluated through summation over
   // quadrature points.
   //
   // We initialy sum all contributions of the
   // residual in the positive sense, so that
   // we don't need to negative the Jacobian
   // entries.  Then, when we sum into the
   // <code>right_hand_side</code> vector,
   // we negate this residual.
   for (unsigned int i=0; i<dofs_per_cell; ++i)
   {
      Sacado::Fad::DFad<double> F_i = 0;
      
      const unsigned int
      component_i = fe_v.get_fe().system_to_component_index(i).first;
      
      // The residual for each row (i) will be accumulating
      // into this fad variable.  At the end of the assembly
      // for this row, we will query for the sensitivities
      // to this variable and add them into the Jacobian.
      
      for (unsigned int point=0; point<n_q_points; ++point)
      {
         if (parameters.is_stationary == false)
            F_i += dWdt[point][component_i] *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);
         
         for (unsigned int d=0; d<dim; d++)
            F_i -= flux[point][component_i][d] *
                   fe_v.shape_grad_component(i, point, component_i)[d] *
                   fe_v.JxW(point);
         
         // Diffusion term for shocks
         if(parameters.diffusion_coef > 0.0)
            for (unsigned int d=0; d<dim; d++)
               F_i += mu_shock(cell_no) *
                      grad_W[point][component_i][d] *
                      fe_v.shape_grad_component(i, point, component_i)[d] *
                      fe_v.JxW(point);
         
         F_i -= parameters.gravity *
                forcing[point][component_i] *
                fe_v.shape_value_component(i, point, component_i) *
                fe_v.JxW(point);
      }
      
      // At the end of the loop, we have to
      // add the sensitivities to the
      // matrix and subtract the residual
      // from the right hand side. Trilinos
      // FAD data type gives us access to
      // the derivatives using
      // <code>F_i.fastAccessDx(k)</code>,
      // so we store the data in a
      // temporary array. This information
      // about the whole row of local dofs
      // is then added to the Trilinos
      // matrix at once (which supports the
      // data types we have chosen).
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
         local_matrix (i, k) += F_i.fastAccessDx(k);
      }
      local_vector (i) -= F_i.val();
   }
   
   delete[] forcing;
   delete[] flux;
   
}


//------------------------------------------------------------------------------
// Contribution from boundary faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_boundary_term (DoFInfo& dinfo, 
                                                    CellInfo& info)
{
   FullMatrix<double>& local_matrix = dinfo.matrix(0).matrix;
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;
   const unsigned int& face_no = dinfo.face_number;
   const double& face_diameter = dinfo.face->diameter();
   const unsigned int& boundary_id = dinfo.face->boundary_id();
   
   const FEValuesBase<dim>& fe_v = info.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   std::vector<Sacado::Fad::DFad<double> >
   independent_local_dof_values (dofs_per_cell);
   
   const unsigned int n_independent_variables = dofs_per_cell;
   
   // Get cell number of two cells adjacent to this face
   const unsigned int cell_no = cell_number (dinfo.cell);
   
   for (unsigned int i = 0; i < dofs_per_cell; i++)
   {
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
      independent_local_dof_values[i].diff(i, n_independent_variables);
   }
   
   // Conservative variable value at face
   Table<2,Sacado::Fad::DFad<double> >
   Wplus (n_q_points, EulerEquations<dim>::n_components),
   Wminus (n_q_points, EulerEquations<dim>::n_components);

   /* TODO:ADIFF
   Table<3,Sacado::Fad::DFad<double> >
   grad_Wplus  (n_q_points, EulerEquations<dim>::n_components, dim);
   */
   
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         Sacado::Fad::DFad<double> dof_theta = parameters.theta *
                                               independent_local_dof_values[i]
                                             +
                                               (1-parameters.theta) *
                                               old_solution(dof_indices[i]);
         Wplus[q][c] += dof_theta *
                        fe_v.shape_value_component(i, q, c);

         /* TODO:ADIFF
         for (unsigned int d = 0; d < dim; d++)
            grad_Wplus[q][c][d] += dof_theta *
                                   fe_v.shape_grad_component(i, q, c)[d];
         */
      }
   

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

   for (unsigned int q = 0; q < n_q_points; q++)
      EulerEquations<dim>::compute_Wminus (boundary_kind,
                                           fe_v.normal_vector(q),
                                           Wplus[q],
                                           boundary_values[q],
                                           Wminus[q]);
   
   
   // Compute numerical flux at all quadrature points
   typedef Sacado::Fad::DFad<double> NormalFlux[EulerEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q], 
                            Wminus[q],
                            cell_average[cell_no],
                            cell_average[cell_no],
                            normal_fluxes[q]);
   
   // Now assemble the face term
   for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true) // comment this for TODO:ADIFF
      {
         Sacado::Fad::DFad<double> F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += normal_fluxes[point][component_i] *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);

            /* TODO:ADIFF
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
         
         for (unsigned int k=0; k<dofs_per_cell; ++k)
            local_matrix (i,k) += F_i.fastAccessDx(k);
         
         local_vector (i) -= F_i.val();
      }
   
   delete[] normal_fluxes;   

}



//------------------------------------------------------------------------------
// Contribution from interior faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                                                CellInfo& info1, CellInfo& info2)
{
   FullMatrix<double>& local_matrix11 = dinfo1.matrix(0,false).matrix;
   FullMatrix<double>& local_matrix12 = dinfo1.matrix(0,true).matrix;
   Vector<double>& local_vector   = dinfo1.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo1.indices;
   const unsigned int& face_no = dinfo1.face_number;
   const double& face_diameter = dinfo1.face->diameter();
   
   FullMatrix<double>& local_matrix_neighbor22 = dinfo2.matrix(0,false).matrix;
   FullMatrix<double>& local_matrix_neighbor21 = dinfo2.matrix(0,true).matrix;
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

   double mu_shock_avg = 0.5 * 
      (mu_shock(cell_no)/dinfo1.cell->diameter()
       +
       mu_shock(neighbor_cell_no)/dinfo2.cell->diameter());

   std::vector<Sacado::Fad::DFad<double> >
   independent_local_dof_values (dofs_per_cell),
   independent_neighbor_dof_values (dofs_per_cell_neighbor);
   
   const unsigned int n_independent_variables = dofs_per_cell + 
                                                dofs_per_cell_neighbor;
   
   for (unsigned int i = 0; i < dofs_per_cell; i++)
   {
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
      independent_local_dof_values[i].diff(i, n_independent_variables);
   }
   
   for (unsigned int i = 0; i < dofs_per_cell_neighbor; i++)
   {
      independent_neighbor_dof_values[i] = current_solution(dof_indices_neighbor[i]);
      independent_neighbor_dof_values[i].diff(i+dofs_per_cell, n_independent_variables);
   }
   
   
   // Compute two states on the face
   Table<2,Sacado::Fad::DFad<double> >
   Wplus (n_q_points, EulerEquations<dim>::n_components),
   Wminus (n_q_points, EulerEquations<dim>::n_components);

   /* TODO:ADIFF
   Table<3,Sacado::Fad::DFad<double> >
   grad_Wplus  (n_q_points, EulerEquations<dim>::n_components, dim),
   grad_Wminus (n_q_points, EulerEquations<dim>::n_components, dim);
   */

   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         Sacado::Fad::DFad<double> dof_theta = parameters.theta *
                                               independent_local_dof_values[i]
                                             +
                                               (1-parameters.theta) *
                                               old_solution(dof_indices[i]);
         Wplus[q][c] += dof_theta * fe_v.shape_value_component(i, q, c);

         /* TODO:ADIFF
         for (unsigned int d = 0; d < dim; d++)
            grad_Wplus[q][c][d] += dof_theta *
                                   fe_v.shape_grad_component(i, q, c)[d];
         */
      }
   
   // Wminus is Neighbouring cell value
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      {
         const unsigned int c = fe_v_neighbor.get_fe().system_to_component_index(i).first;
         Sacado::Fad::DFad<double> dof_theta = parameters.theta *
                                               independent_neighbor_dof_values[i]
                                             +
                                               (1-parameters.theta) *
                                               old_solution(dof_indices_neighbor[i]);
         Wminus[q][c] += dof_theta *
                         fe_v_neighbor.shape_value_component(i, q, c);

         /* TODO:ADIFF
         for (unsigned int d = 0; d < dim; d++)
            grad_Wminus[q][c][d] += dof_theta *
                                    fe_v_neighbor.shape_grad_component(i, q, c)[d];
         */
      }
   
   
   // Compute numerical flux at all quadrature points
   typedef Sacado::Fad::DFad<double> NormalFlux[EulerEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q],
                            Wminus[q],
                            cell_average[cell_no],
                            cell_average[neighbor_cell_no],
                            normal_fluxes[q]);
   
   // Now assemble the face term
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true) // comment for TODO:ADIFF
      {
         Sacado::Fad::DFad<double> F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += normal_fluxes[point][component_i] *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);

            /* TODO:ADIFF
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
         
         for (unsigned int k=0; k<dofs_per_cell; ++k)
            local_matrix11 (i,k) += F_i.fastAccessDx(k);

            for (unsigned int k=0; k<dofs_per_cell_neighbor; ++k)
               local_matrix12(i,k) += F_i.fastAccessDx(dofs_per_cell+k);
         
         local_vector(i) -= F_i.val();
      }
   
   // Contributions to neighbouring cell
   for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      if (fe_v_neighbor.get_fe().has_support_on_face(i, face_no_neighbor) == true) // comment for TODO:ADIFF
      {
         Sacado::Fad::DFad<double> F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            
            F_i -= normal_fluxes[point][component_i] *
                   fe_v_neighbor.shape_value_component(i, point, component_i) *
                   fe_v_neighbor.JxW(point);

            /* TODO:ADIFF
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
         
         for (unsigned int k=0; k<dofs_per_cell_neighbor; ++k)
            local_matrix_neighbor22 (i,k) += F_i.fastAccessDx(dofs_per_cell+k);
         
         for (unsigned int k=0; k<dofs_per_cell; ++k)
            local_matrix_neighbor21 (i,k) += F_i.fastAccessDx(k);
         
         local_vector_neighbor (i) -= F_i.val();
      }
   
   delete[] normal_fluxes;

}

//------------------------------------------------------------------------------
// Assemble matrices
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::assemble_system (IntegratorImplicit<dim>& integrator)
{
   // Compute artificial viscosity for shock capturing
   compute_mu_shock ();

   system_matrix = 0;
   right_hand_side = 0;

   MeshWorker::loop<dim,dim,MeshWorker::DoFInfo<dim>,MeshWorker::IntegrationInfoBox<dim> >
   (dof_handler.begin_active(),
    dof_handler.end(),
    integrator.dof_info, 
    integrator.info_box,
    boost::bind(&ConservationLaw<dim>::integrate_cell_term, 
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::integrate_boundary_term,
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::integrate_face_term,
                this, _1, _2, _3, _4),
    integrator.assembler);
}

template class ConservationLaw<2>;
