#include <base/quadrature_lib.h>

#include <fe/fe_values.h>
#include <fe/fe_values_extractors.h>

#include <dofs/dof_handler.h>

#include "equation.h"
#include "claw.h"

using namespace dealii;

//-----------------------------------------------------------------------------
// Positivity limiter of Zhang-Shu
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_positivity_limiter ()
{
   TimerOutput::Scope t(computing_timer, "Positivity limiter");

   if(fe.degree == 0) return;
   
   const double gas_gamma = EulerEquations<dim>::gas_gamma;
   const unsigned int density_component = EulerEquations<dim>::density_component;
   const unsigned int energy_component = EulerEquations<dim>::energy_component;
   
   // Find mininimum density and pressure in the whole grid
   const double eps_tol = 1.0e-13;
   double eps = eps_tol;
   {
      typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active(),
         endc = dof_handler.end();
      for (; cell!=endc; ++cell)
      if(cell->is_locally_owned())
      {
         const unsigned int c = cell_number (cell);

         eps = std::min(eps, cell_average[c][density_component]);
         double pressure = EulerEquations<dim>::template compute_pressure<double> (cell_average[c]);
         eps = std::min(eps, pressure);
         if(eps < eps_tol)
         {
            std::cout << "\n Negative state at position " << cell->center() << "\n\n";
            AssertThrow(false, ExcMessage("Fatal: Negative states"));
         }
      }
      eps = -Utilities::MPI::max(-eps, mpi_communicator);
   }
   
   // Need 2N - 3 >= degree for the quadrature to be exact.
   // Choose same order as used for assembly process.
   unsigned int N = fe.degree + 2;
   QGaussLobatto<dim>  quadrature_formula(N);
   const unsigned int n_q_points = quadrature_formula.size();
   FEValues<dim> fe_values (mapping(), fe, quadrature_formula, update_values);
   std::vector<double> density_values(n_q_points), energy_values(n_q_points);
   std::vector< Tensor<1,dim> > momentum_values(n_q_points);
   std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
   
   const FEValuesExtractors::Scalar density  (density_component);
   const FEValuesExtractors::Scalar energy   (energy_component);
   const FEValuesExtractors::Vector momentum (0);
   
   right_hand_side = current_solution;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for(; cell != endc; ++cell)
   if(cell->is_locally_owned())
   {
      unsigned int c = cell_number(cell);
      fe_values.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      
      // First limit density
      fe_values[density].get_function_values(right_hand_side, density_values);
      
      // find minimum density at GLL points
      double rho_min = 1.0e20;
      for(unsigned int q=0; q<n_q_points; ++q)
         rho_min = std::min(rho_min, density_values[q]);
      
      double density_average = cell_average[c][density_component];
      double rat = std::fabs(density_average - eps) /				//is it eps or eps_tol?
                   (std::fabs(density_average - rho_min) + 1.0e-13);
      double theta1 = std::min(rat, 1.0);
      
      if(theta1 < 1.0)
      {
         if(parameters.basis == Parameters::AllParameters<dim>::Qk)
         {
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int comp_i = fe.system_to_component_index(i).first;
               if(comp_i == density_component)
                  right_hand_side(local_dof_indices[i]) =
                     theta1         * right_hand_side(local_dof_indices[i])
                     + (1.0 - theta1) * density_average;
            }
         }
         else
         {
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int comp_i = fe.system_to_component_index(i).first;
               unsigned int base_i = fe.system_to_component_index(i).second;
               if(comp_i == density_component && base_i > 0)
                  right_hand_side(local_dof_indices[i]) *= theta1;
            }
         }
      }
      
      // now limit pressure
      fe_values[density].get_function_values(right_hand_side, density_values);
      fe_values[momentum].get_function_values(right_hand_side, momentum_values);
      fe_values[energy].get_function_values(right_hand_side, energy_values);
      
      double energy_average = cell_average[c][energy_component];
      Tensor<1,dim> momentum_average;
      for(unsigned int i=0; i<dim; ++i)
         momentum_average[i] = cell_average[c][i];
      
      double theta2 = 1.0;
      for(unsigned int q=0; q<n_q_points; ++q)
      {
         double pressure = (gas_gamma-1.0)*(energy_values[q] -
                                            0.5*momentum_values[q].norm_square()/density_values[q]);
         if(pressure < eps)
         {
            double drho = density_values[q] - density_average;
            Tensor<1,dim> dm = momentum_values[q] - momentum_average;
            double dE = energy_values[q] - energy_average;
            double a1 = 2.0*drho*dE - dm*dm;
            double b1 = 2.0*drho*(energy_average - eps/(gas_gamma-1.0))
                      + 2.0*density_average*dE
                      - 2.0*momentum_average*dm;
            double c1 = 2.0*density_average*energy_average
                      - momentum_average*momentum_average
                      - 2.0*eps*density_average/(gas_gamma-1.0);
            // Divide by a1 to avoid round-off error
            b1 /= a1; c1 /= a1;
            double D = std::sqrt( std::fabs(b1*b1 - 4.0*c1) );
            double t1 = 0.5*(-b1 - D);
            double t2 = 0.5*(-b1 + D);
            double t;
            if(t1 > -1.0e-12 && t1 < 1.0 + 1.0e-12)
               t = t1;
            else if(t2 > -1.0e-12 && t2 < 1.0 + 1.0e-12)
               t = t2;
            else
            {
               std::cout << "Problem in positivity limiter\n";
               std::cout << "\t a1, b1, c1 = " << a1 << " " << b1 << " " << c1 << "\n";
               std::cout << "\t t1, t2 = " << t1 << " " << t2 << "\n";
               std::cout << "\t eps, rho_min = " << eps << " " << rho_min << "\n";
               std::cout << "\t theta1 = " << theta1 << "\n";
               std::cout << "\t pressure = " << pressure << "\n";
               exit(0);
            }
            // t should strictly lie in [0,1]
            t = std::min(1.0, t);
            t = std::max(0.0, t);
            // Need t < 1.0. If t==1 upto machine precision
            // then we are suffering from round off error.
            // In this case we take the cell average value, t=0.
            if(std::fabs(1.0-t) < 1.0e-14) t = 0.0;
            theta2 = std::min(theta2, t);
         }
      }
      
      if(theta2 < 1.0)
      {
         if(parameters.basis == Parameters::AllParameters<dim>::Qk)
         {
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int comp_i = fe.system_to_component_index(i).first;
               right_hand_side(local_dof_indices[i]) =
                  theta2         * right_hand_side(local_dof_indices[i])
                  + (1.0 - theta2)  * cell_average[c][comp_i];
            }
         }
         else
         {
            for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
               unsigned int base_i = fe.system_to_component_index(i).second;
               if(base_i > 0)
                  right_hand_side(local_dof_indices[i]) *= theta2;
            }
            
         }
      }
      
   }
   current_solution = right_hand_side;
}

//-----------------------------------------------------------------------------
// Positivity limiter of Zhang-Shu for one cell
// This function is called from within TVD limiter.
// Limit the solution in newton_update.
//-----------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::apply_positivity_limiter_cell
   (typename DoFHandler<dim>::active_cell_iterator& cell,
    PosLimData<dim>& data)
{   
   static const double gas_gamma = EulerEquations<dim>::gas_gamma;
   static const unsigned int density_component = EulerEquations<dim>::density_component;
   static const unsigned int energy_component  = EulerEquations<dim>::energy_component;
   
   // Find mininimum density and pressure in the whole grid
   static const double eps = 1.0e-13;
   
   // Need 2N - 3 >= degree for the quadrature to be exact.
   // Choose same order as used for assembly process.
   QGaussLobatto<dim>&  quadrature_formula = data.quadrature_formula;
   const unsigned int& n_q_points = data.n_q_points;
   FEValues<dim>& fe_values = data.fe_values;
   std::vector<double>& density_values = data.density_values;
   std::vector<double>& energy_values = data.energy_values;
   std::vector< Tensor<1,dim> >& momentum_values = data.momentum_values;
   std::vector<unsigned int>& local_dof_indices = data.local_dof_indices;
   
   static const FEValuesExtractors::Scalar density (density_component);
   static const FEValuesExtractors::Scalar energy  (energy_component);
   static const FEValuesExtractors::Vector momentum(0);
   
   unsigned int c = cell_number(cell);
   fe_values.reinit(cell);
   
   // First limit density
   fe_values[density].get_function_values(newton_update, density_values);
   
   // find minimum density at GLL points
   double rho_min = 1.0e20;
   for(unsigned int q=0; q<n_q_points; ++q)
      rho_min = std::min(rho_min, density_values[q]);
   
   double density_average = cell_average[c][density_component];
   double rat = std::fabs(density_average - eps) /
               (std::fabs(density_average - rho_min) + 1.0e-13);
   double theta1 = std::min(rat, 1.0);
   
   if(theta1 < 1.0)
   {
      cell->get_dof_indices (local_dof_indices);
      if(parameters.basis == Parameters::AllParameters<dim>::Qk)
      {
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         {
            unsigned int comp_i = fe.system_to_component_index(i).first;
            if(comp_i == density_component)
               newton_update(local_dof_indices[i]) =
                  theta1           * newton_update(local_dof_indices[i])
                  + (1.0 - theta1) * density_average;
         }
      }
      else
      {
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         {
            unsigned int comp_i = fe.system_to_component_index(i).first;
            unsigned int base_i = fe.system_to_component_index(i).second;
            if(comp_i == density_component && base_i > 0)
               newton_update(local_dof_indices[i]) *= theta1;
         }
      }
   }
   
   // now limit pressure
   fe_values[density].get_function_values(newton_update, density_values);
   fe_values[momentum].get_function_values(newton_update, momentum_values);
   fe_values[energy].get_function_values(newton_update, energy_values);
   
   double energy_average = cell_average[c][energy_component];
   Tensor<1,dim> momentum_average;
   for(unsigned int i=0; i<dim; ++i)
      momentum_average[i] = cell_average[c][i];
   
   double theta2 = 1.0;
   for(unsigned int q=0; q<n_q_points; ++q)
   {
      double pressure = (gas_gamma-1.0)*(energy_values[q] -
                                         0.5*momentum_values[q].norm_square()/density_values[q]);
      if(pressure < eps)
      {
         double drho = density_values[q] - density_average;
         Tensor<1,dim> dm = momentum_values[q] - momentum_average;
         double dE = energy_values[q] - energy_average;
         double a1 = 2.0*drho*dE - dm*dm;
         double b1 = 2.0*drho*(energy_average - eps/(gas_gamma-1.0))
                     + 2.0*density_average*dE
                     - 2.0*momentum_average*dm;
         double c1 = 2.0*density_average*energy_average
                     - momentum_average*momentum_average
                     - 2.0*eps*density_average/(gas_gamma-1.0);
         // Divide by a1 to avoid round-off error
         b1 /= a1; c1 /= a1;
         double D = std::sqrt( std::fabs(b1*b1 - 4.0*c1) );
         double t1 = 0.5*(-b1 - D);
         double t2 = 0.5*(-b1 + D);
         double t;
         if(t1 > -1.0e-12 && t1 < 1.0 + 1.0e-12)
            t = t1;
         else if(t2 > -1.0e-12 && t2 < 1.0 + 1.0e-12)
            t = t2;
         else
         {
            std::cout << "Problem in positivity limiter\n";
            std::cout << "\t a1, b1, c1 = " << a1 << " " << b1 << " " << c1 << "\n";
            std::cout << "\t t1, t2 = " << t1 << " " << t2 << "\n";
            std::cout << "\t eps, rho_min = " << eps << " " << rho_min << "\n";
            std::cout << "\t theta1 = " << theta1 << "\n";
            std::cout << "\t pressure = " << pressure << "\n";
            exit(0);
         }
         // t should strictly lie in [0,1]
         t = std::min(1.0, t);
         t = std::max(0.0, t);
         // Need t < 1.0. If t==1 upto machine precision
         // then we are suffering from round off error.
         // In this case we take the cell average value, t=0.
         if(std::fabs(1.0-t) < 1.0e-14) t = 0.0;
         theta2 = std::min(theta2, t);
      }
   }
   
   if(theta2 < 1.0)
   {
      if(!(theta1<1.0)) // local_dof_indices has not been computed before
         cell->get_dof_indices (local_dof_indices);
      if(parameters.basis == Parameters::AllParameters<dim>::Qk)
      {
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         {
            unsigned int comp_i = fe.system_to_component_index(i).first;
            newton_update(local_dof_indices[i]) =
               theta2            * newton_update(local_dof_indices[i])
               + (1.0 - theta2)  * cell_average[c][comp_i];
         }
      }
      else
      {
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         {
            unsigned int base_i = fe.system_to_component_index(i).second;
            if(base_i > 0)
               newton_update(local_dof_indices[i]) *= theta2;
         }
         
      }
   }
}

template class ConservationLaw<2>;
