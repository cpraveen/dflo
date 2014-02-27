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
      case Parameters::Limiter::grad:
         apply_limiter_grad ();
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
   
   QTrapez<1>        q_trapez;
   QMidpoint<1>      q_midpoint;
   QAnisotropic<dim> qrule_x (q_trapez, q_midpoint);
   FEValues<dim>     fe_values_x (fe, qrule_x, update_values);
   QAnisotropic<dim> qrule_y (q_midpoint, q_trapez);
   FEValues<dim>     fe_values_y (fe, qrule_y, update_values);
   
   Quadrature<dim> qsupport (fe.get_unit_support_points());
   FEValues<dim>   fe_values (fe, qsupport, update_q_points);
   
   std::vector<Vector<double> > face_values_x(2,
                                              Vector<double>(EulerEquations<dim>::n_components));
   std::vector<Vector<double> > face_values_y(2,
                                              Vector<double>(EulerEquations<dim>::n_components));
    
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
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell0 = dof_handler0.begin_active(),
      endc0 = dof_handler0.end();
   
   for(; cell != endc; ++cell, ++cell0)
   {
      const unsigned int c = cell_number(cell);
      const double dx = cell->diameter() / std::sqrt(2.0);
      const double Mdx2 = parameters.M * dx * dx;
      
      cell0->get_dof_indices (cell_indices);
      
      // Backward difference of cell averages
      dbx = 0;
      if(lcell[c] != endc0)
      {
         get_cell_average (lcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            dbx(i) = cell_average(cell_indices[i]) - avg_nbr(i);
      }
      
      // Forward difference of cell averages
      dfx = 0;
      if(rcell[c] != endc0)
      {
         get_cell_average (rcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            dfx(i) = avg_nbr(i) - cell_average(cell_indices[i]);
      }
       
      // Backward difference of cell averages
      dby = 0;
      if(bcell[c] != endc0)
      {
         get_cell_average (bcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            dby(i) = cell_average(cell_indices[i]) - avg_nbr(i);
      }
      
      // Forward difference of cell averages
      dfy = 0;
      if(tcell[c] != endc0)
      {
         get_cell_average (tcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            dfy(i) = avg_nbr(i) - cell_average(cell_indices[i]);
      }
      
      fe_values_x.reinit(cell);
      fe_values_x.get_function_values(current_solution, face_values_x);
      fe_values_y.reinit(cell);
      fe_values_y.get_function_values(current_solution, face_values_y);
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         Dx(i) = face_values_x[1][i] - face_values_x[0][i];
         Dy(i) = face_values_y[1][i] - face_values_y[0][i];
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
            unsigned int base_i = fe.system_to_component_index(i).second;
            Point<dim> dr = p[base_i] - cell->center();
            current_solution(dof_indices[i]) = cell_average(cell_indices[comp_i])
               + dr[0] * Dx_new(comp_i) + dr[1] * Dy_new(comp_i);
         }
      }
      
   }
}

//------------------------------------------------------------------------------
// Apply gradient limiter
// Note: This is implemented only for 2-D.
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_limiter_grad ()
{
   if(fe.degree == 0) return;
   const unsigned int n_components = EulerEquations<dim>::n_components;

   QGauss<dim> qrule (fe.degree + 1);
   FEValues<dim> fe_values_grad (fe, qrule, update_gradients | update_JxW_values);

   Quadrature<dim> qsupport (fe.get_unit_support_points());
   FEValues<dim>   fe_values (fe, qsupport, update_q_points);
   
   Vector<double> dfx (n_components);
   Vector<double> dbx (n_components);
   Vector<double> Dx  (n_components);
   
   Vector<double> dfy (n_components);
   Vector<double> dby (n_components);
   Vector<double> Dy  (n_components);
   
   Vector<double> Dx_new (n_components);
   Vector<double> Dy_new (n_components);
   Vector<double> avg_nbr (n_components);
   
   std::vector<unsigned int> cell_indices (fe0.dofs_per_cell);
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   std::vector< std::vector< Tensor<1,dim> > > grad (qrule.size(), std::vector< Tensor<1,dim> >(n_components));
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell0 = dof_handler0.begin_active(),
      endc0 = dof_handler0.end();
   
   for(; cell != endc; ++cell, ++cell0)
   {
      const unsigned int c = cell_number(cell);
      const double dx = cell->diameter() / std::sqrt(2.0);
      const double Mdx2 = parameters.M * dx * dx;
      
      cell0->get_dof_indices (cell_indices);
      
      // Backward difference of cell averages
      dbx = 0;
      if(lcell[c] != endc0)
      {
         get_cell_average (lcell[c], avg_nbr);
         for(unsigned int i=0; i<n_components; ++i)
            dbx(i) = cell_average(cell_indices[i]) - avg_nbr(i);
      }
      
      // Forward difference of cell averages
      dfx = 0;
      if(rcell[c] != endc0)
      {
         get_cell_average (rcell[c], avg_nbr);
         for(unsigned int i=0; i<n_components; ++i)
            dfx(i) = avg_nbr(i) - cell_average(cell_indices[i]);
      }
      
      // Backward difference of cell averages
      dby = 0;
      if(bcell[c] != endc0)
      {
         get_cell_average (bcell[c], avg_nbr);
         for(unsigned int i=0; i<n_components; ++i)
            dby(i) = cell_average(cell_indices[i]) - avg_nbr(i);
      }
      
      // Forward difference of cell averages
      dfy = 0;
      if(tcell[c] != endc0)
      {
         get_cell_average (tcell[c], avg_nbr);
         for(unsigned int i=0; i<n_components; ++i)
            dfy(i) = avg_nbr(i) - cell_average(cell_indices[i]);
      }
      
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
         Dx(i) = 0.5 * dx * avg_grad[0];
         Dy(i) = 0.5 * dx * avg_grad[1];
      }
      
      // Transform to characteristic variables
      typedef double EigMatrix[n_components][n_components];
      EigMatrix Rx, Lx, Ry, Ly;
      if(parameters.char_lim)
      {
         Vector<double> avg (n_components);
         for(unsigned int i=0; i<n_components; ++i)
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
      for(unsigned int i=0; i<n_components; ++i)
      {
         Dx_new(i) = minmod(Dx(i), dbx(i), dfx(i), Mdx2);
         Dy_new(i) = minmod(Dy(i), dby(i), dfy(i), Mdx2);
         change_x += std::fabs(Dx_new(i) - Dx(i));
         change_y += std::fabs(Dy_new(i) - Dy(i));
      }
      change_x /= n_components;
      change_y /= n_components;
      
      // If limiter is active, reduce polynomial to linear
      if(change_x + change_y > 1.0e-10)
      {
         Dx_new /= 0.5 * dx;
         Dy_new /= 0.5 * dx;
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
            unsigned int base_i = fe.system_to_component_index(i).second;
            Point<dim> dr = p[base_i] - cell->center();
            current_solution(dof_indices[i]) = cell_average(cell_indices[comp_i])
               + dr[0] * Dx_new(comp_i) + dr[1] * Dy_new(comp_i);
         }
      }
      
   }
}

template class ConservationLaw<2>;
