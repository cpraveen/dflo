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
   TimerOutput::Scope t(computing_timer, "Limiter");

   if(parameters.basis == Parameters::AllParameters<dim>::Qk)
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
   else
   {
      switch(parameters.limiter_type)
      {
         case Parameters::Limiter::none:
            break;
         case Parameters::Limiter::TVB:
            apply_limiter_TVB_Pk ();
            break;
         default:
            AssertThrow(false, ExcMessage("Unknown limiter_type"));
      }
      
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
   
   // Data for positivity limiter
   PosLimData<dim> pos_lim_data (fe, mapping());
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      endc0 = dh_cell.end();
   
   const double beta = parameters.beta;
   newton_update = current_solution;
   
   for(; cell != endc; ++cell)
   if(cell->is_locally_owned())
   {
      const unsigned int c = cell_number(cell);
      if(shock_indicator[c] > 1.0)
      {
         const double dx = cell->diameter() / std::sqrt(1.0*dim);
         const double Mdx2 = parameters.M * dx * dx;
         
         // Compute average gradient in cell
         fe_values_grad.reinit(cell);
         fe_values_grad.get_function_gradients(newton_update, grad);
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
         if(lcell[c] != endc0)
         {
            get_cell_average (lcell[c], avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dbx(i) = cell_average[c][i] - avg_nbr(i);
         }
         
         // Forward difference of cell averages
         dfx = Dx;
         if(rcell[c] != endc0)
         {
            get_cell_average (rcell[c], avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dfx(i) = avg_nbr(i) - cell_average[c][i];
         }
         
         // Backward difference of cell averages
         dby = Dy;
         if(bcell[c] != endc0)
         {
            get_cell_average (bcell[c], avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dby(i) = cell_average[c][i] - avg_nbr(i);
         }
         
         // Forward difference of cell averages
         dfy = Dy;
         if(tcell[c] != endc0)
         {
            get_cell_average (tcell[c], avg_nbr);
            for(unsigned int i=0; i<n_components; ++i)
               dfy(i) = avg_nbr(i) - cell_average[c][i];
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
               newton_update(dof_indices[i]) = cell_average[c][comp_i]
                                                  + dr[0] * Dx_new(comp_i)
                                                  + dr[1] * Dy_new(comp_i);
            }
         }
      }
      if(parameters.pos_lim)
         apply_positivity_limiter_cell (cell, pos_lim_data);
   }
   current_solution = newton_update;
}

//------------------------------------------------------------------------------
// Apply TVB limiter
// Note: This is implemented only for 2-D.
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_limiter_TVB_Pk ()
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
   
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   
   static const double sqrt_3 = sqrt(3.0);
   const double beta = 0.5 * parameters.beta;

   // Data for positivity limiter
   PosLimData<dim> pos_lim_data (fe, mapping());

   newton_update = current_solution;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      endc0 = dh_cell.end();
   
   for(; cell != endc; ++cell)
   if(cell->is_locally_owned())
   {
      const unsigned int c = cell_number(cell);
      if(shock_indicator[c] > 1.0)
      {
         const double dx = cell->diameter() / std::sqrt(1.0*dim);
         const double Mdx2 = parameters.M * dx * dx;
         
         cell->get_dof_indices(dof_indices);
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         {
            unsigned int comp_i = fe.system_to_component_index(i).first;
            unsigned int base_i = fe.system_to_component_index(i).second;
            if(base_i == 1)
               Dx(comp_i) = newton_update(dof_indices[i]) * sqrt_3;
            else if(base_i == fe.degree+1)
               Dy(comp_i) = newton_update(dof_indices[i]) * sqrt_3;
         }
         
         // angular momentum for square cells = v_x - u_y
         const double ang_mom = Dx(1) - Dy(0);
         
         // Backward difference of cell averages
         dbx = Dx;
         if(lcell[c] != endc0)
         {
            get_cell_average (lcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dbx(i) = cell_average[c][i] - avg_nbr(i);
         }
         
         // Forward difference of cell averages
         dfx = Dx;
         if(rcell[c] != endc0)
         {
            get_cell_average (rcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dfx(i) = avg_nbr(i) - cell_average[c][i];
         }
         
         // Backward difference of cell averages
         dby = Dy;
         if(bcell[c] != endc0)
         {
            get_cell_average (bcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dby(i) = cell_average[c][i] - avg_nbr(i);
         }
         
         // Forward difference of cell averages
         dfy = Dy;
         if(tcell[c] != endc0)
         {
            get_cell_average (tcell[c], avg_nbr);
            for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
               dfy(i) = avg_nbr(i) - cell_average[c][i];
         }
         
         // Transform to characteristic variables
         typedef double EigMatrix[EulerEquations<dim>::n_components][EulerEquations<dim>::n_components];
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
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
         {
            Dx_new(i) = minmod(Dx(i), beta*dbx(i), beta*dfx(i), Mdx2);
            Dy_new(i) = minmod(Dy(i), beta*dby(i), beta*dfy(i), Mdx2);
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
            if(parameters.conserve_angular_momentum)
            {
               Dy_new(0) = 0.5 * (Dy_new(0) - (ang_mom - Dx_new(1)));
               Dx_new(1) = ang_mom + Dy_new(0);
            }
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int comp_i = fe.system_to_component_index(i).first;
               unsigned int base_i = fe.system_to_component_index(i).second;
               if(base_i == 1)
                  newton_update(dof_indices[i]) = Dx_new(comp_i) / sqrt_3;
               else if(base_i == fe.degree + 1)
                  newton_update(dof_indices[i]) = Dy_new(comp_i) / sqrt_3;
               else if(base_i != 0)
                  newton_update(dof_indices[i]) = 0.0;
            }
         }
         
      }
      if(parameters.pos_lim)
         apply_positivity_limiter_cell (cell, pos_lim_data);
   }

   current_solution = newton_update;
}

template class ConservationLaw<2>;
