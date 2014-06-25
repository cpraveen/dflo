#include <base/quadrature_lib.h>

#include <fe/fe_values.h>

#include <dofs/dof_handler.h>

#include "equation.h"
#include "claw.h"

using namespace dealii;

//------------------------------------------------------------------------------
// TVB version of minmod limiter. If Mdx2=0 then it is TVD limiter.
//------------------------------------------------------------------------------
double minmod (const double& a,
               const double& b,
               const double& c,
               const double& Mdx2)
{
   double aa = std::fabs(a);
   if(aa < Mdx2) return a;
   
   if(a*b > 0 && b*c > 0)
   {
      double s = (a > 0) ? 1.0 : -1.0;
      return s * std::min(aa, std::min(std::fabs(b), std::fabs(c)));
   }
   else
      return 0;
}

//------------------------------------------------------------------------------
// Apply selected limiter
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_limiter ()
{
   switch(parameters.limiter_type)
   {
      case Parameters::Limiter::none:
         break;
      case Parameters::Limiter::TVB:
         apply_limiter_TVB ();
         break;
   }
}

//------------------------------------------------------------------------------
// Apply TVB limiter
// Note: This is implemented only for 2-D.
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_limiter_TVB ()
{
   if(fe.degree == 0) return;
   
   Vector<double> dfx (EulerEquations<dim>::n_components);
   Vector<double> dbx (EulerEquations<dim>::n_components);
   Vector<double> Dx  (EulerEquations<dim>::n_components);
    
   Vector<double> dfy (EulerEquations<dim>::n_components);
   Vector<double> dby (EulerEquations<dim>::n_components);
   Vector<double> Dy  (EulerEquations<dim>::n_components);
    
   Vector<double> Dx_new (EulerEquations<dim>::n_components);
   Vector<double> Dy_new (EulerEquations<dim>::n_components);
   Vector<double> avg_nbr (EulerEquations<dim>::n_components);
   
   std::vector<unsigned int> cell_indices (fe0.dofs_per_cell);
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);

   static const double sqrt_3 = sqrt(3.0);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell0 = dof_handler0.begin_active(),
      endc0 = dof_handler0.end();
   
   for(; cell != endc; ++cell, ++cell0)
   {
      const unsigned int c = cell_number(cell);
      if(shock_indicator[c] > 1)
      {
         const double dx = cell->diameter() / std::sqrt(2.0);
         const double Mdx2 = parameters.M * dx * dx;
         
         cell0->get_dof_indices (cell_indices);
         
         cell->get_dof_indices(dof_indices);
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         {
            unsigned int comp_i = fe.system_to_component_index(i).first;
            unsigned int base_i = fe.system_to_component_index(i).second;
            if(base_i == 1)
               Dx(comp_i) = current_solution(dof_indices[i]) * sqrt_3;
            else if(base_i == fe.degree+1)
               Dy(comp_i) = current_solution(dof_indices[i]) * sqrt_3;
         }
         
         // Backward difference of cell averages
         dbx = Dx;
         if(lcell[c] != endc0)
         {
            get_cell_average (lcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dbx(i) = cell_average(cell_indices[i]) - avg_nbr(i);
         }
         
         // Forward difference of cell averages
         dfx = Dx;
         if(rcell[c] != endc0)
         {
            get_cell_average (rcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dfx(i) = avg_nbr(i) - cell_average(cell_indices[i]);
         }
         
         // Backward difference of cell averages
         dby = Dy;
         if(bcell[c] != endc0)
         {
            get_cell_average (bcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dby(i) = cell_average(cell_indices[i]) - avg_nbr(i);
         }
         
         // Forward difference of cell averages
         dfy = Dy;
         if(tcell[c] != endc0)
         {
            get_cell_average (tcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dfy(i) = avg_nbr(i) - cell_average(cell_indices[i]);
         }
         
         // Transform to characteristic variables
         typedef double EigMatrix[EulerEquations<dim>::n_components][EulerEquations<dim>::n_components];
         EigMatrix Rx, Lx, Ry, Ly;
         if(parameters.char_lim)
         {
            Vector<double> avg (EulerEquations<dim>::n_components);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               avg(i) = cell_average(cell_indices[i]);
            EulerEquations<dim>::compute_eigen_matrix (avg, Rx, Lx, Ry, Ly);
            EulerEquations<dim>::transform_to_char (Lx, dbx);
            EulerEquations<dim>::transform_to_char (Lx, dfx);
            EulerEquations<dim>::transform_to_char (Ly, dby);
            EulerEquations<dim>::transform_to_char (Ly, dfy);
            EulerEquations<dim>::transform_to_char (Lx, Dx);
            EulerEquations<dim>::transform_to_char (Ly, Dy);
         }
         
         // Apply minmod limiter
         double change_x = 0;
         double change_y = 0;
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
         {
            Dx_new(i) = minmod(Dx(i), dbx(i), dfx(i), Mdx2);
            Dy_new(i) = minmod(Dy(i), dby(i), dfy(i), Mdx2);
            change_x += std::fabs(Dx_new(i) - Dx(i));
            change_y += std::fabs(Dy_new(i) - Dy(i));
         }
         change_x /= EulerEquations<dim>::n_components;
         change_y /= EulerEquations<dim>::n_components;
         
         // If limiter is active, reduce polynomial to linear
         if(change_x + change_y > 1.0e-10)
         {
            if(parameters.char_lim)
            {
               EulerEquations<dim>::transform_to_con (Rx, Dx_new);
               EulerEquations<dim>::transform_to_con (Ry, Dy_new);
            }
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int comp_i = fe.system_to_component_index(i).first;
               unsigned int base_i = fe.system_to_component_index(i).second;
               if(base_i == 1)
                  current_solution(dof_indices[i]) = Dx_new(comp_i) / sqrt_3;
               else if(base_i == fe.degree + 1)
                  current_solution(dof_indices[i]) = Dy_new(comp_i) / sqrt_3;
               else if(base_i != 0)
                  current_solution(dof_indices[i]) = 0.0;
            }
         }
         
      }
   }
}

template class ConservationLaw<2>;
