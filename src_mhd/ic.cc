#include "ic.h"
#include "equation.h"
#include "claw.h"

using namespace dealii;

//--------------------------------------------------------------------------------------------
// Initial condition for Rayleigh-Taylor problem
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void RayleighTaylor<dim>::vector_value (const Point<dim> &p,
                                        Vector<double>   &values) const
{
   // Density
   if(p[1] < 0.0)
      values[MHDEquations<dim>::density_component] = 1.0;
   else
      values[MHDEquations<dim>::density_component] = 2.0;
   
   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[d] = 0.0;
   
   double vel = A *
                (1.0 + std::cos(2.0*numbers::PI*p[0]/Lx))/2.0 *
                (1.0 + std::cos(2.0*numbers::PI*p[1]/Ly))/2.0;
   
   values[1] = values[MHDEquations<dim>::density_component] * vel;
   
   double pressure = P0 - gravity * values[MHDEquations<dim>::density_component] * p[1];

   // Energy
   values[MHDEquations<dim>::energy_component] =
      pressure/(MHDEquations<dim>::gas_gamma - 1.0)
      + 0.5 * values[MHDEquations<dim>::density_component] * vel * vel;
}

//--------------------------------------------------------------------------------------------
// Initial condition for isentropic vortex problem
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void IsentropicVortex<dim>::vector_value (const Point<dim> &p,
                                          Vector<double>   &values) const	
{
   const double gamma = MHDEquations<dim>::gas_gamma;
   
   double r2   = (p[0]-x0)*(p[0]-x0) + (p[1]-y0)*(p[1]-y0);
   
   double rho = std::pow(1.0 - a2*exp(1.0-r2), 1.0/(gamma-1.0));
   double vex =  - a1 * (p[1]-y0) * std::exp(0.5*(1.0-r2));
   double vey =  + a1 * (p[0]-x0) * std::exp(0.5*(1.0-r2));
   double pre = std::pow(rho, gamma) / gamma;
   
   values[0] = rho * vex;
   values[1] = rho * vey;
   values[MHDEquations<dim>::density_component] = rho;
   values[MHDEquations<dim>::energy_component] = pre/(gamma-1.0)
                                                   + 0.5 * rho * (vex*vex + vey*vey);
}

//--------------------------------------------------------------------------------------------
// Initial condition for system of three isentropic vortices.
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void VortexSystem<dim>::vector_value (const Point<dim> &p,
                                      Vector<double>   &values) const
{
   const double gamma = MHDEquations<dim>::gas_gamma;
   
   double rho = 0, vex = 0, vey = 0, pre = 0;
   for(unsigned int i=0; i<3; ++i)
   {
      double r2   = (p[0]-x[i])*(p[0]-x[i]) + (p[1]-y[i])*(p[1]-y[i]);
      
      rho += std::pow(1.0 - a2*exp(1.0-r2), 1.0/(gamma-1.0));
      vex +=  - a1 * (p[1]-y[i]) * std::exp(0.5*(1.0-r2));
      vey +=  + a1 * (p[0]-x[i]) * std::exp(0.5*(1.0-r2));
   }
   
   rho -= 2.0;
   vex /= 3.0;
   vey /= 3.0;
   pre = std::pow(rho, gamma) / gamma;
   
   // Put large pressure in the region { |x| < 0.1 and |y| < 0.1 }
   if(std::fabs(p[0]) < 0.1 && std::fabs(p[1]) < 0.1) pre = 50.0;
   
   values[MHDEquations<dim>::momentum_component] = rho * vex;
   values[MHDEquations<dim>::momentum_component+1] = rho * vey;
   values[MHDEquations<dim>::density_component] = rho;
   values[MHDEquations<dim>::energy_component] = pre/(gamma-1.0)
                                                 + 0.5 * rho * (vex*vex + vey*vey);
   for(unsigned int i=0; i<dim; ++i)
     values[MHDEquations<dim>::magnetic_component+i]=0;
}

//--------------------------------------------------------------------------------------------
// Keplerian disk
// TODO TO BE COMPLETED
//--------------------------------------------------------------------------------------------
template <int dim>
void KeplerianDisk<dim>::vector_value (const Point<dim> &p,
                                       Vector<double>   &values) const
{
   const double gamma = MHDEquations<dim>::gas_gamma;
   
   double r   = p.norm();
   double vtheta = 1/std::sqrt(r);
   
   double vex = -vtheta * p[1] / r;
   double vey =  vtheta * p[0] / r;
   
   double rho = 1.0;

   if(r < r0-rs || r > r1+rs)
      rho = rho_out;
   else if(r >= r0-rs && r <= r1+rs)
      rho = rho_disk;
   else if(r <= r0+rs)
      rho = 0;
   else if(r >= r1-rs)
      rho = 0;
   
   values[0] = rho * vex;
   values[1] = rho * vey;
   values[MHDEquations<dim>::density_component] = rho;
   values[MHDEquations<dim>::energy_component] = pressure/(gamma-1.0)
                                                   + 0.5 * rho * (vex*vex + vey*vey);
}

//------------------------------------------------------------------------------
// Sets initial condition based on input file.
// For Qk basis we can just do interpolation.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_Qk ()
{
   if(parameters.ic_function == "rt")
      VectorTools::interpolate(mapping(), dof_handler,
                               RayleighTaylor<dim>(parameters.gravity), old_solution);
   else if(parameters.ic_function == "isenvort")
      VectorTools::interpolate(mapping(), dof_handler,
                               IsentropicVortex<dim>(5.0, 0.0, 0.0), old_solution);
   else if(parameters.ic_function == "vortsys")
      VectorTools::interpolate(mapping(), dof_handler,
                               VortexSystem<dim>(), old_solution);
   else
      VectorTools::interpolate(mapping(), dof_handler,
                               parameters.initial_conditions, old_solution);
   
   for(unsigned int i = 0; i < old_solution.size(); i++)
     if(isnan(old_solution[i]))
       std::cout<<"\n \t The initial caondition has a NaN values! in : "<< i << "\n" ;
   
   current_solution = old_solution;
   predictor = old_solution;
}

//------------------------------------------------------------------------------
// Sets initial condition based on input file.
// For Pk basis we have to do an L2 projection.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_Pk ()
{
   RayleighTaylor<dim> rayleigh_taylor(parameters.gravity);
   IsentropicVortex<dim> isentropic_vortex(5.0, 0.0, 0.0);
   VortexSystem<dim> vortex_system;
   Function<dim>* ic_function;
   if(parameters.ic_function == "rt")
      ic_function = &rayleigh_taylor;
   else if(parameters.ic_function == "isenvort")
      ic_function = &isentropic_vortex;
   else if(parameters.ic_function == "vortsys")
      ic_function = &vortex_system;
   else
      ic_function = &parameters.initial_conditions;

   QGauss<dim> quadrature (fe.degree+1);
   unsigned int n_q_points = quadrature.size();
   FEValues<dim> fe_values (mapping(), fe, quadrature, 
                            update_values|update_q_points|update_JxW_values);
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   
   unsigned int n_comp=4;
   
   /*if(parameters.equation=='euler')
     n_comp=EulerEquations<dim>::n_components;
   if(parameters.equation=='mhd')
     n_comp=MHDEquations<dim>::n_components;//*/

   std::vector< Vector<double> > ic_values(n_q_points, Vector<double>(n_comp));
   
   old_solution = 0.0;

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      fe_values.reinit (cell);
      ic_function->vector_value_list(fe_values.get_quadrature_points(), ic_values);
      cell->get_dof_indices(dof_indices);
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
         unsigned int comp_i = fe.system_to_component_index(i).first;
         for(unsigned int q=0; q<n_q_points; ++q)
            old_solution(dof_indices[i]) += ic_values[q][comp_i] *
                                            fe_values.shape_value_component(i, q, comp_i) *
                                            fe_values.JxW(q);
         }

      unsigned int c = cell_number(cell);
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         old_solution(dof_indices[i]) *= inv_mass_matrix[c][i];
   }
   
   current_solution = old_solution;
   predictor = old_solution;

}

//------------------------------------------------------------------------------
// Set intitial condition by interpolation or projection depending on
// type of basis function.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition ()
{
   pcout << "Setting initial condition\n";
   
   if(parameters.basis == Parameters::AllParameters<dim>::Qk)
      set_initial_condition_Qk();
   else
      set_initial_condition_Pk();
}

template class RayleighTaylor<2>;
template class IsentropicVortex<2>;
template class VortexSystem<2>;
template class ConservationLaw<2>;
