#include <base/quadrature_lib.h>

#include <fe/fe_values.h>

#include <dofs/dof_handler.h>

#include "equation.h"
#include "claw.h"

using namespace dealii;

//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_limiter_TVB ()
{
   if(fe.degree == 0) return;
   
   QTrapez<1> q_trapez;
   QMidpoint<1> q_midpoint;
   QAnisotropic<dim> qrule_x (q_trapez, q_midpoint);
   FEValues<dim> fe_values_x (fe, qrule_x, update_values);
   QAnisotropic<dim> qrule_y (q_midpoint, q_trapez);
   FEValues<dim> fe_values_y (fe, qrule_y, update_values);
   
   Quadrature<dim> qsupport (fe.get_unit_support_points());
   FEValues<dim> fe_values (fe, qsupport, update_q_points);
   
   std::vector<Vector<double> > face_values(2,
                                            Vector<double>(EulerEquations<dim>::n_components));
   Vector<double> df (EulerEquations<dim>::n_components);
   Vector<double> db (EulerEquations<dim>::n_components);
   Vector<double> DF (EulerEquations<dim>::n_components);
   Vector<double> DB (EulerEquations<dim>::n_components);
   Vector<double> DFx_new (EulerEquations<dim>::n_components);
   Vector<double> DBx_new (EulerEquations<dim>::n_components);
   Vector<double> DFy_new (EulerEquations<dim>::n_components);
   Vector<double> DBy_new (EulerEquations<dim>::n_components);
   Vector<double> avg_nbr (EulerEquations<dim>::n_components);
   
   std::vector<unsigned int> cell_indices (fe0.dofs_per_cell);
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell0 = dof_handler0.begin_active(),
      endc0 = dof_handler0.end();
   
   for(unsigned int c=0; cell != endc; ++c, ++cell, ++cell0)
   {
      double dx = cell->diameter() / std::sqrt(2.0);
      double Mdx2 = parameters.M * dx * dx;
      
      cell0->get_dof_indices (cell_indices);
      
      // Limit x derivative
      fe_values_x.reinit(cell);
      fe_values_x.get_function_values(current_solution, face_values);

      // Backward difference of cell averages
      db = 0;
      if(lcell[c] != endc0)
      {
         get_cell_average (lcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            db(i) = cell_average(cell_indices[i]) - avg_nbr(i);
      }
      
      // Forward difference of cell averages
      df = 0;
      if(rcell[c] != endc0)
      {
         get_cell_average (rcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            df(i) = avg_nbr(i) - cell_average(cell_indices[i]);
      }
      
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         DB(i) = cell_average(cell_indices[i]) - face_values[0][i];
         DF(i) = face_values[1][i] - cell_average(cell_indices[i]);
      }
      
      // Transform to characteristic variables
      if(parameters.char_lim)
      {
         
      }
      
      // Apply minmod limiter
      double change_x = 0;
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         DBx_new(i) = minmod(DB(i), db(i), df(i), Mdx2);
         DFx_new(i) = minmod(DF(i), db(i), df(i), Mdx2);
         change_x += std::fabs(DBx_new(i) - DB(i)) + std::fabs(DFx_new(i) - DF(i));
      }
      change_x /= 2 * EulerEquations<dim>::n_components;
      
      // Limit y derivative
      fe_values_y.reinit(cell);
      fe_values_y.get_function_values(current_solution, face_values);
      
      // Backward difference of cell averages
      db = 0;
      if(bcell[c] != endc0)
      {
         get_cell_average (bcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            db(i) = cell_average(cell_indices[i]) - avg_nbr(i);
      }
      
      // Forward difference of cell averages
      df = 0;
      if(tcell[c] != endc0)
      {
         get_cell_average (tcell[c], avg_nbr);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
            df(i) = avg_nbr(i) - cell_average(cell_indices[i]);
      }
      
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         DB(i) = cell_average(cell_indices[i]) - face_values[0][i];
         DF(i) = face_values[1][i] - cell_average(cell_indices[i]);
      }
      
      // Transform to characteristic variables
      if(parameters.char_lim)
      {
         
      }
      
      // Apply minmod limiter
      double change_y = 0;
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         DBy_new(i) = minmod(DB(i), db(i), df(i), Mdx2);
         DFy_new(i) = minmod(DF(i), db(i), df(i), Mdx2);
         change_y += std::fabs(DBy_new(i) - DB(i)) + std::fabs(DFy_new(i) - DF(i));
      }
      change_y /= 2 * EulerEquations<dim>::n_components;
      
      // If limiter is active, reduce polynomial to linear
      if(change_x + change_y > 1.0e-10)
      {
         Vector<double> gradx(EulerEquations<dim>::n_components);
         Vector<double> grady(EulerEquations<dim>::n_components);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
         {
            gradx(i) = 0.5*(DBx_new(i) + DFx_new(i))/dx;
            grady(i) = 0.5*(DBy_new(i) + DFy_new(i))/dx;
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
               + dr[0] * gradx(comp_i) + dr[1] * grady(comp_i);
         }
      }
      
   }
}

template class ConservationLaw<2>;