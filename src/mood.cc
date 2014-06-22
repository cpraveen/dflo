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

#include <base/quadrature_lib.h>
#include <fe/fe_values.h>
#include <dofs/dof_handler.h>

#include "equation.h"
#include "claw.h"

using namespace dealii;

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

//--------------------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::reduce_degree(const typename DoFHandler<dim>::cell_iterator &cell,
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
         else
         {
            re_update_new[c] = true;
            
            if(cell_degree[c] > 0)
            {
               // restrict solution
               reduce_degree(cell, c, fe_values);
               terminate = false;
            }
            else
            {
               unsigned int c = 0;
               // Loop over neighbours of cell
               for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
               {
                  if (cell->at_boundary(f) == false)
                     if (cell->neighbor(f)->has_children() == false || cell->neighbor_is_coarser(f))
                     {
                        unsigned int cn = cell_number(cell->neighbor(f));
                        if(cell_degree[cn] > 0 && re_update_new[cn]==false)
                        {
                           // restrict solution
                           reduce_degree(cell->neighbor(f), cn, fe_values);
                           re_update_new[cn] = true;
                           terminate = false;
                           ++c;
                        }
                     }
                     else
                     {
                        for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
                        {
                           unsigned int cn = cell_number(cell->neighbor_child_on_subface (f, subface));
                           if(cell_degree[cn] > 0 && re_update_new[cn]==false)
                           {
                              // restrict solution
                              reduce_degree(cell->neighbor_child_on_subface(f,subface), cn, fe_values);
                              re_update_new[cn] = true;
                              terminate = false;
                              ++c;
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
