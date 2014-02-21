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
   Vector<double> DFx (EulerEquations<dim>::n_components);
   Vector<double> DBx (EulerEquations<dim>::n_components);
    
   Vector<double> dfy (EulerEquations<dim>::n_components);
   Vector<double> dby (EulerEquations<dim>::n_components);
   Vector<double> DFy (EulerEquations<dim>::n_components);
   Vector<double> DBy (EulerEquations<dim>::n_components);
    
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
   
   for(; cell != endc; ++cell, ++cell0)
   {
      unsigned int c = cell_number(cell);
      double dx = cell->diameter() / std::sqrt(2.0);
      double Mdx2 = parameters.M * dx * dx;
      
      cell0->get_dof_indices (cell_indices);
      
      // Limit x derivative
      fe_values_x.reinit(cell);
      fe_values_x.get_function_values(current_solution, face_values_x);

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
       
      // Limit y derivative
      fe_values_y.reinit(cell);
      fe_values_y.get_function_values(current_solution, face_values_y);
      
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
      
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         DBx(i) = cell_average(cell_indices[i]) - face_values_x[0][i];
         DFx(i) = face_values_x[1][i] - cell_average(cell_indices[i]);
         DBy(i) = cell_average(cell_indices[i]) - face_values_y[0][i];
         DFy(i) = face_values_y[1][i] - cell_average(cell_indices[i]);
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
         EulerEquations<dim>::transform_to_char (Lx, DBx);
         EulerEquations<dim>::transform_to_char (Lx, DFx);
         EulerEquations<dim>::transform_to_char (Ly, DBy);
         EulerEquations<dim>::transform_to_char (Ly, DFy);
      }
      
      // Apply minmod limiter
      double change_x = 0;
      double change_y = 0;
      for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
      {
         DBx_new(i) = minmod(DBx(i), dbx(i), dfx(i), Mdx2);
         DFx_new(i) = minmod(DFx(i), dbx(i), dfx(i), Mdx2);
         DBy_new(i) = minmod(DBy(i), dby(i), dfy(i), Mdx2);
         DFy_new(i) = minmod(DFy(i), dby(i), dfy(i), Mdx2);
         change_x += std::fabs(DBx_new(i) - DBx(i)) + std::fabs(DFx_new(i) - DFx(i));
         change_y += std::fabs(DBy_new(i) - DBy(i)) + std::fabs(DFy_new(i) - DFy(i));
      }
      change_x /= 2 * EulerEquations<dim>::n_components;
      change_y /= 2 * EulerEquations<dim>::n_components;
      
      // If limiter is active, reduce polynomial to linear
      if(change_x + change_y > 1.0e-10)
      {
         Vector<double> gradx(EulerEquations<dim>::n_components);
         Vector<double> grady(EulerEquations<dim>::n_components);
         for(unsigned int i=0; i<EulerEquations<dim>::n_components; ++i)
         {
            gradx(i) = (DBx_new(i) + DFx_new(i))/dx;
            grady(i) = (DBy_new(i) + DFy_new(i))/dx;
         }
         if(parameters.char_lim)
         {
            EulerEquations<dim>::transform_to_con (Rx, gradx);
            EulerEquations<dim>::transform_to_con (Ry, grady);
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
