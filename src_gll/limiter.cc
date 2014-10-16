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
            apply_limiter_TVB_Qk ();
            break;
         default:
            AssertThrow(false, ExcMessage("Unknown limiter_type"));
      }
 }

//------------------------------------------------------------------------------
// Apply gradient limiter
// Note: This is implemented only for 2-D.
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_limiter_TVB_Qk ()
{
   if(fe.degree == 0) return;
   const unsigned int n_components = EulerEquations<dim>::n_components;

   QGauss<dim> qrule (fe.degree + 1);
   FEValues<dim> fe_values_grad (mapping(), fe, qrule, update_gradients | update_JxW_values);

   // NOTE: We get multiple sets of same support points since fe is an FESystem
   Quadrature<dim> qsupport (fe.get_unit_support_points());
   FEValues<dim>   fe_values (mapping(), fe, qsupport, update_q_points);
   
   Vector<double> dfx (n_components);
   Vector<double> dbx (n_components);
   Vector<double> Dx  (n_components);
   
   Vector<double> dfy (n_components);
   Vector<double> dby (n_components);
   Vector<double> Dy  (n_components);
   
   Vector<double> Dx_new (n_components);
   Vector<double> Dy_new (n_components);
   Vector<double> avg_nbr (n_components);
   
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   std::vector< std::vector< Tensor<1,dim> > > grad (qrule.size(),
                                                     std::vector< Tensor<1,dim> >(n_components));
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      endc0 = dh_cell.end();
   
   const double beta = parameters.beta;
   
   for(; cell != endc; ++cell)
   {
      const unsigned int c = cell_number(cell);
      if(shock_indicator[c] > 1.0)
      {
         const double dx = cell->diameter() / std::sqrt(1.0*dim);
         const double Mdx2 = parameters.M * dx * dx;
         
         // Compute average gradient in cell
         fe_values_grad.reinit(cell);
         fe_values_grad.get_function_gradients(current_solution, grad);
         Tensor<1,dim> avg_grad;
         
         for(unsigned int i=0; i<n_components; ++i)
         {
            avg_grad = 0;
            for(unsigned int q=0; q<qrule.size(); ++q)
               avg_grad += grad[q][i] * fe_values_grad.JxW(q);
            avg_grad /= cell->measure();
            Dx(i) = dx * avg_grad[0];
            Dy(i) = dx * avg_grad[1];
         }
         
         // Backward difference of cell averages
         dbx = Dx;
         if(cell_data[c].lcell != endc0)
         {
            get_cell_average (cell_data[c].lcell, avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dbx(i) = cell_average[c][i] - avg_nbr(i);
         }
         else if(cell_data[c].lbc == EulerEquations<dim>::no_penetration_boundary)
         {
            dbx    = 0.0;
            dbx(0) = 2.0*cell_average[c][0];
         }
         
         // Forward difference of cell averages
         dfx = Dx;
         if(cell_data[c].rcell != endc0)
         {
            get_cell_average (cell_data[c].rcell, avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dfx(i) = avg_nbr(i) - cell_average[c][i];
         }
         else if(cell_data[c].rbc == EulerEquations<dim>::no_penetration_boundary)
         {
            dfx    = 0.0;
            dfx(0) = -2.0*cell_average[c][0];
         }
         
         // Backward difference of cell averages
         dby = Dy;
         if(cell_data[c].bcell != endc0)
         {
            get_cell_average (cell_data[c].bcell, avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dby(i) = cell_average[c][i] - avg_nbr(i);
         }
         else if(cell_data[c].bbc == EulerEquations<dim>::no_penetration_boundary)
         {
            dby    = 0.0;
            dby(1) = 2.0*cell_average[c][1];
         }
         
         // Forward difference of cell averages
         dfy = Dy;
         if(cell_data[c].tcell != endc0)
         {
            get_cell_average (cell_data[c].tcell, avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dfy(i) = avg_nbr(i) - cell_average[c][i];
         }
         else if(cell_data[c].tbc == EulerEquations<dim>::no_penetration_boundary)
         {
            dfy    = 0.0;
            dfy(1) = -2.0*cell_average[c][1];
         }
         
         // Transform to characteristic variables
         typedef double EigMatrix[n_components][n_components];
         EigMatrix Rx, Lx, Ry, Ly;
         if(parameters.char_lim)
         {
            EulerEquations<dim>::compute_eigen_matrix (cell_average[c], Rx, Lx, Ry, Ly);
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
         for(unsigned int i=0; i<n_components; ++i)
         {
            Dx_new(i) = minmod(Dx(i), beta*dbx(i), beta*dfx(i), Mdx2);
            Dy_new(i) = minmod(Dy(i), beta*dby(i), beta*dfy(i), Mdx2);
            change_x += std::fabs(Dx_new(i) - Dx(i));
            change_y += std::fabs(Dy_new(i) - Dy(i));
         }
         change_x /= n_components;
         change_y /= n_components;
         
         // If limiter is active, reduce polynomial to linear
         if(change_x + change_y > 1.0e-10)
         {
            Dx_new /= dx;
            Dy_new /= dx;
            if(parameters.char_lim)
            {
               EulerEquations<dim>::transform_to_con (Rx, Dx_new);
               EulerEquations<dim>::transform_to_con (Ry, Dy_new);
            }
            cell->get_dof_indices(dof_indices);
            fe_values.reinit (cell);
            const std::vector<Point<dim> >& p = fe_values.get_quadrature_points();
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int comp_i = fe.system_to_component_index(i).first;
               Point<dim> dr = p[i] - cell->center();
               current_solution(dof_indices[i]) = cell_average[c][comp_i]
                                                  + dr[0] * Dx_new(comp_i)
                                                  + dr[1] * Dy_new(comp_i);
            }
         }
      }
   }
}

template class ConservationLaw<2>;
