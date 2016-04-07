/*
for rk=1,2,3,...
   find min/max
   re_update(1:ncell) = true
   cell_degree(1:ncell) = degree
   terminate = false
   while(! terminate)
      assemble cells with re_update==true
      update cells with re_update==true
      terminate = true
      for cell=1,2,... and with re_update(cell)==true
         check DMP
         if DMP is satisfied 
            set re_update(cell) = false
         else
            if(cell_degree(cell) > 0)
               cell_degree(cell) -= 1
               restrict solution to cell_degree(cell)
               terminate = false
            else 
               for each neighbor of cell
                  if(cell_degree(neighbor) > 0 and re_update(neighbor)==false)
                     cell_degree(neighbor) -= 1
                     re_update(neighbor) = true
                     restrict solution to cell_degree(neighbor)
                     terminate = false
                  endif
               endfor
            endif
         endif
   endwhile
   apply positivity limiter
endfor
*/

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>

#include "equation.h"
#include "claw.h"

using namespace dealii;

//--------------------------------------------------------------------------------------------
// For Qk basis, compute degree reduction matrices
//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::compute_reduction_matrices()
{
   if(fe.degree == 0) return;
   
   std::cout << "Constructing interpolation matrices ...\n";
   const FE_DGQArbitraryNodes<dim> forig(QGauss<1>(fe.degree+1));
   unsigned int ndofs_orig = forig.dofs_per_cell;
   
   Rmatrix.resize(fe.degree, FullMatrix<double>(ndofs_orig,ndofs_orig));
   
   for(int deg=0; deg<fe.degree; ++deg)
   {
      const FE_DGQArbitraryNodes<dim> ff(QGauss<1>(deg+1));
      unsigned int ndofs = ff.dofs_per_cell;
      FullMatrix<double> m1(ndofs, ndofs_orig);
      FullMatrix<double> m2(ndofs_orig, ndofs);
      ff.get_interpolation_matrix(forig, m1);
      forig.get_interpolation_matrix(ff, m2);
      // Rmatrix[deg] = m2*m1
      m2.mmult(Rmatrix[deg], m1);
   }
}

//--------------------------------------------------------------------------------------------
// For each cell, compute min and max density in {cell, neighbors}
//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::compute_min_max_mood_var()
{
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for(; cell != endc; ++cell)
   {
      unsigned int c = cell_number(cell);
      min_mood_var[c] = cell_average[c][EulerEquations<dim>::density_component];
      max_mood_var[c] = cell_average[c][EulerEquations<dim>::density_component];

      // Loop over neighbours of cell
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
         if (cell->at_boundary(f) == false)
         {
            if (cell->neighbor(f)->has_children() == false || cell->neighbor_is_coarser(f))
            {
               unsigned int cn = cell_number(cell->neighbor(f));
               min_mood_var[cn] = std::min( min_mood_var[cn],
                                            cell_average[cn][EulerEquations<dim>::density_component]);
               max_mood_var[cn] = std::max( max_mood_var[cn],
                                            cell_average[cn][EulerEquations<dim>::density_component]);
            }
            else
            {
               for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
               {
                  unsigned int cn = cell_number(cell->neighbor_child_on_subface (f, subface));
                  min_mood_var[cn] = std::min( min_mood_var[cn],
                                               cell_average[cn][EulerEquations<dim>::density_component]);
                  max_mood_var[cn] = std::max( max_mood_var[cn],
                                               cell_average[cn][EulerEquations<dim>::density_component]);
               }
            }
         }
      }
   }
}

//--------------------------------------------------------------------------------------------
// Reduce degree of Qk solution. Interpolate to lower degree and interpolate back
// to original degree. Adjust to preserve cell average value.
//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::reduce_degree_Qk(const typename DoFHandler<dim>::cell_iterator &cell,
                                            const unsigned int cell_no,
                                            FEValues<dim> &fe_values)
{
   --cell_degree[cell_no];
   
   Vector<double> reduced_solution(fe.dofs_per_cell);
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   Vector<double> cell_avg_reduced(EulerEquations<dim>::n_components);
   
   cell->get_dof_indices (dof_indices);
   reduced_solution = 0;
   cell_avg_reduced = 0;
   
   fe_values.reinit(cell);
   
   // no of dofs in each component of fe. we just take from component 0 since
   // all components are same for us.
   unsigned int ndofs = fe.base_element(0).dofs_per_cell;
   const unsigned int deg = cell_degree[cell_no];
   
   for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
   {
      unsigned int comp_i = fe.system_to_component_index(i).first;
      unsigned int base_i = fe.system_to_component_index(i).second;
      for(unsigned int j=0; j<ndofs; ++j)
      {
         unsigned int index_j = fe.component_to_system_index(comp_i,j);
         reduced_solution(i) += Rmatrix[deg](base_i,j) * current_solution(dof_indices[index_j]);
      }
      
      cell_avg_reduced(comp_i) += reduced_solution(i) * fe_values.JxW(base_i);
   }
   
   // compute cell average of reduced_solution
   cell_avg_reduced /= cell->measure();
   
   // adjust cell average
   for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
   {
      unsigned int comp_i = fe.system_to_component_index(i).first;
      current_solution(dof_indices[i]) = cell_average[cell_no][comp_i]
                                       + reduced_solution(i)
                                       - cell_avg_reduced(comp_i);
   }
   
}

//--------------------------------------------------------------------------------------------
// Reduce degree of Pk solution by one. All coefficients higher than this degree are set
// to zero.
//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::reduce_degree_Pk(const typename DoFHandler<dim>::cell_iterator &cell,
                                            const unsigned int cell_no,
                                            FEValues<dim> &fe_values)
{
   --cell_degree[cell_no];
   
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   cell->get_dof_indices (dof_indices);
   for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
   {
      unsigned int base_i = fe.system_to_component_index(i).second;
      unsigned int deg = index_to_degree[base_i];
      if(deg > cell_degree[cell_no])
         current_solution(dof_indices[i]) = 0.0;
   }
   
}

//--------------------------------------------------------------------------------------------
// Reduce degree of solution by one
//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::reduce_degree(const typename DoFHandler<dim>::cell_iterator &cell,
                                         const unsigned int cell_no,
                                         FEValues<dim> &fe_values)
{
   if(parameters.basis == Parameters::AllParameters<dim>::Qk)
      reduce_degree_Qk(cell, cell_no, fe_values);
   else
      reduce_degree_Pk(cell, cell_no, fe_values);
}

//--------------------------------------------------------------------------------------------
// Returns second order coefficients of density in its legendre expansion
// NOTE: We should divide by dx*dx etc. for adaptive meshes.
//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::get_mood_second_derivatives
   (const typename DoFHandler<dim>::cell_iterator &cell,
    std::vector<double>& D2)
{
   Assert(dim == 2, ExcNotImplemented());
   Assert(parameters.basis == Parameters::AllParameters<dim>::Pk,
          ExcNotImplemented());
   
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   cell->get_dof_indices(dof_indices);
   
   // This gets rho_xx and rho_yy
   unsigned int ndofs = fe.dofs_per_cell / EulerEquations<dim>::n_components;
   unsigned int shift = ndofs * EulerEquations<dim>::density_component;
   D2[0] = current_solution(dof_indices[shift+2]);
   D2[1] = current_solution(dof_indices[shift+2*fe.degree+1]);
}

//--------------------------------------------------------------------------------------------
// Apply u2 test of Diot on cell. If satisfied, return true, else return false.
//--------------------------------------------------------------------------------------------
template <int dim>
bool ConservationLaw<dim>::test_u2(const typename DoFHandler<dim>::cell_iterator &cell)
{
   // If degree < 2, then we should not apply this test.
   // Hence always return false.
   if(fe.degree < 2) return false;
   
   std::vector<double> D2(dim);
   get_mood_second_derivatives (cell, D2);
   
   std::vector<double> D2_min(D2), D2_max(D2);
   
   for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
   {
      if (cell->at_boundary(f) == false)
      {
         if (cell->neighbor(f)->has_children() == false || cell->neighbor_is_coarser(f))
         {
            get_mood_second_derivatives (cell->neighbor(f), D2);
            for(unsigned int d=0; d<dim; ++d)
            {
               D2_min[d] = std::min(D2_min[d], D2[d]);
               D2_max[d] = std::max(D2_max[d], D2[d]);
            }
         }
         else
         {
            for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
            {
               get_mood_second_derivatives (cell->neighbor_child_on_subface(f,subface), D2);
               for(unsigned int d=0; d<dim; ++d)
               {
                  D2_min[d] = std::min(D2_min[d], D2[d]);
                  D2_max[d] = std::max(D2_max[d], D2[d]);
               }
            }
         }
      }
   }
   
   static const double eps = 0.5;
   static const double fact = 1.0 - eps;
   for(unsigned int d=0; d<dim; ++d)
   {
      if(D2_min[d]*D2_max[d] < 0) return false;
      if(std::fabs(D2_min[d]) < std::fabs(D2_max[d]) * fact) return false;
   }
   
   return true;
}
//--------------------------------------------------------------------------------------------
// For each cell, if re_update = true, then check DMP property. If not satisfied, reduce
// polynomial degree.
//
// If DMP property is not satisfied for any cell and its degree = 0, then we reduce the
// degree of all its neighbouring cells.
//--------------------------------------------------------------------------------------------
template <int dim>
bool ConservationLaw<dim>::apply_mood(unsigned int &n_reduce,
                                      unsigned int &n_re_update,
                                      unsigned int &n_reset)
{
   if(fe.degree == 0) return true;
   
   QGauss<dim> quadrature(fe.degree+1);
   FEValues<dim> fe_values(mapping(), fe, quadrature, update_JxW_values);
   const double eps = 1.0e-6;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   std::vector<bool> re_update_new (triangulation.n_active_cells());
   
   bool terminate = true;
   for(; cell != endc; ++cell)
   {
      unsigned int c = cell_number(cell);
      
      if(re_update[c] == true && re_update_new[c] == false)
      {
         // check DMP property
         // min_mood_var <= mood_var <= max_mood_var
         const double mood_var = cell_average[c][EulerEquations<dim>::density_component];
         const bool dmp = (mood_var - min_mood_var[c] > -eps && mood_var - max_mood_var[c] < eps);
         
         if(dmp)
         {
            re_update_new[c] = false;
         }
         else if(test_u2(cell) == true)
         {
            re_update_new[c] = false;
         }
         else
         {
            if(cell_degree[c] > 1)
            {
               // restrict solution
               reduce_degree(cell, c, fe_values);
               re_update_new[c] = true;
               terminate = false;
            }
            else if(cell_degree[c] == 1 && shock_indicator[c] < 1.0)
            {
               shock_indicator[c] = 1.0e20; // switch on limiter for this cell
               re_update_new[c] = false;
            }
            else
            {
               unsigned int c = 0;
               // Loop over neighbours of cell
               for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
               {
                  if (cell->at_boundary(f) == false)
                  {
                     if (cell->neighbor(f)->has_children() == false || cell->neighbor_is_coarser(f))
                     {
                        unsigned int cn = cell_number(cell->neighbor(f));
                        if(cell_degree[cn] > 1 && re_update_new[cn]==false)
                        {
                           // restrict solution
                           reduce_degree(cell->neighbor(f), cn, fe_values);
                           re_update_new[cn] = true;
                           terminate = false;
                           ++c;
                        }
                        else if(cell_degree[cn] == 1 &&
                                re_update_new[cn]==false &&
                                shock_indicator[cn] < 1.0)
                        {
                           shock_indicator[cn] = 1.0e20;
                           re_update_new[cn] = false;
                           ++c;
                        }
                     }
                     else
                     {
                        for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
                        {
                           unsigned int cn = cell_number(cell->neighbor_child_on_subface (f, subface));
                           if(cell_degree[cn] > 1 && re_update_new[cn]==false)
                           {
                              // restrict solution
                              reduce_degree(cell->neighbor_child_on_subface(f,subface), cn, fe_values);
                              re_update_new[cn] = true;
                              terminate = false;
                              ++c;
                           }
                           else if(cell_degree[cn] == 1 &&
                                   re_update_new[cn]==false &&
                                   shock_indicator[cn] < 1.0)
                           {
                              shock_indicator[cn] = 1.0e20;
                              re_update_new[cn] = false;
                              ++c;
                           }
                        }
                     }
                  }
               }
               //AssertThrow(c>0, ExcMessage("MOOD failure"));
            }
         }
      }
   }
   
   // Now we update the re_update flag
   for(unsigned int i=0; i<triangulation.n_active_cells(); ++i)
      re_update[i] = false;
   
   n_reduce = 0;

   cell = dof_handler.begin_active();
   for(; cell != endc; ++cell)
   {
      unsigned int c = cell_number (cell);
      
      if(re_update_new[c]) ++n_reduce;
      
      if(re_update[c] == false && re_update_new[c] == true)
      {
         for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
         {
            if (cell->at_boundary(f) == false)
            {
               if (cell->neighbor(f)->has_children() == false || cell->neighbor_is_coarser(f))
               {
                  unsigned int cn = cell_number(cell->neighbor(f));
                  re_update[cn] = true;
               }
               else
               {
                  for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
                  {
                     unsigned int cn = cell_number(cell->neighbor_child_on_subface (f, subface));
                     re_update[cn] = true;
                  }
               }
            }
         }
         re_update[c] = true;
      }
   }
   
   // mark cells whose solution has to be reset to previous stage
   std::vector<bool> reset_sol(triangulation.n_active_cells(), false);
   n_re_update = 0;
   cell = dof_handler.begin_active();
   for(; cell != endc; ++cell)
   {
      unsigned int c = cell_number(cell);
      if(re_update[c])
      {
         // reset solution to previous stage
         ++n_re_update;
         
         reset_sol[c] = true;
         
         for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
         {
            if (cell->at_boundary(f) == false)
            {
               if (cell->neighbor(f)->has_children() == false || cell->neighbor_is_coarser(f))
               {
                  unsigned int cn = cell_number(cell->neighbor(f));
                  reset_sol[cn] = true;
               }
               else
               {
                  for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
                  {
                     unsigned int cn = cell_number(cell->neighbor_child_on_subface (f, subface));
                     reset_sol[cn] = true;
                  }
               }
            }
         }
      }
   }
   
   // now reset the solution to previous stage
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   n_reset = 0;
   cell = dof_handler.begin_active();
   for(; cell != endc; ++cell)
   {
      unsigned int c = cell_number(cell);
      if(reset_sol[c])
      {
         ++n_reset;
         cell->get_dof_indices(dof_indices);
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            current_solution(dof_indices[i]) = predictor(dof_indices[i]);
      }
   }

   return terminate;
}

template class ConservationLaw<2>;
