#include "ic.h"
#include "equation.h"
#include "claw.h"

using namespace dealii;

// Temperatures used in Rayleigh-Taylor problem
const double Tl = 2.0;
const double Tu = 1.0;

//--------------------------------------------------------------------------------------------
template <int dim>
void PolytropicHydrostatic<dim>::vector_value (const Point<dim> &p,
                                               Vector<double>   &values) const
{
   // Density
   values[EulerEquations<dim>::density_component] = std::pow(rho0,nu-1) - alpha * (nu-1.0)/nu * p[1];
   values[EulerEquations<dim>::density_component]
      = std::pow(values[EulerEquations<dim>::density_component], 1.0/(nu-1.0));
   
   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[d] = 0.0;
   
   double pressure = alpha * std::pow(values[EulerEquations<dim>::density_component],nu);
   
   // Energy
   values[EulerEquations<dim>::energy_component] = pressure/(EulerEquations<dim>::gas_gamma - 1.0);
}

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
      values[EulerEquations<dim>::density_component] = 1.0;
   else
      values[EulerEquations<dim>::density_component] = 2.0;
   
   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[d] = 0.0;
   
   double vel = A *
                (1.0 + std::cos(2.0*numbers::PI*p[0]/Lx))/2.0 *
                (1.0 + std::cos(2.0*numbers::PI*p[1]/Ly))/2.0;
   
   values[1] = values[EulerEquations<dim>::density_component] * vel;
   
   double pressure = P0 - gravity * values[EulerEquations<dim>::density_component] * p[1];

   // Energy
   values[EulerEquations<dim>::energy_component] =
      pressure/(EulerEquations<dim>::gas_gamma - 1.0)
      + 0.5 * values[EulerEquations<dim>::density_component] * vel * vel;
}

//--------------------------------------------------------------------------------------------
template <int dim>
void RadialRayleighTaylor<dim>::vector_value (const Point<dim> &p,
                                              Vector<double>   &values) const
{
   double r = p.norm();
   double theta = std::atan2(p[1],p[0]);
   double alpha = std::exp(-r0)/(std::exp(-r0) + drho);
   double factor = std::exp(r0*(1.0-alpha)/alpha);
   
   double pressure;
   
   if(r < r0)
   {
      pressure = std::exp(-r);
   }
   else
   {
      pressure = factor * std::exp(-r/alpha);
   }
   
   double ds = 0.01; // smoothing length = 2*ds
   double ri = r0*(1.0 + eta*cos(k*theta));
   double smoothH = 0.5*(1.0 - std::tanh((r-ri)/ds));
   
   values[EulerEquations<dim>::density_component]
      = std::exp(-r) * smoothH
        + factor/alpha * std::exp(-r/alpha) * (1.0 - smoothH);

   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[d] = 0.0;
   
   // Energy
   values[EulerEquations<dim>::energy_component] =  pressure/(EulerEquations<dim>::gas_gamma - 1.0);
}

//--------------------------------------------------------------------------------------------
// Isothermal hydrostatic test case from Xing and Shu
//--------------------------------------------------------------------------------------------
template <int dim>
void IsothermalHydrostatic<dim>::vector_value (const Point<dim> &p,
                                               Vector<double>   &values) const
{
   double ff1 = - (rho0 * g)/p0 * (p[0] + p[1]);
   double ff2 = - (100 * rho0 * g)/p0 * ((p[0]-0.3)*(p[0]-0.3) + (p[1]-0.3)*(p[1]-0.3));
   
   values[EulerEquations<dim>::density_component] = rho0 * std::exp(ff1);
   
   double pressure = p0 * std::exp(ff1) + eta * std::exp(ff2);
   
   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[d] = 0.0;
   
   // Energy
   values[EulerEquations<dim>::energy_component] =  pressure/(EulerEquations<dim>::gas_gamma - 1.0);
}

//--------------------------------------------------------------------------------------------
// Isothermal hydrostatic test case from Xing and Shu
//--------------------------------------------------------------------------------------------
template <int dim>
void UnsteadyGravity<dim>::vector_value (const Point<dim> &p,
                                               Vector<double>   &values) const
{
   values[EulerEquations<dim>::density_component]
      = 1.0 + 0.2 * std::sin(M_PI*(p[0]+p[1] - time*(u0+v0)));
   
   
   // Momentum
   values[0] = values[EulerEquations<dim>::density_component] * u0;
   values[1] = values[EulerEquations<dim>::density_component] * v0;
   
   double pressure = p0 - p[0] - p[1] + time*(u0+v0)
                     + 0.2 * std::cos(M_PI*(p[0]+p[1]-time*(u0+v0)))/M_PI;

   // Energy
   values[EulerEquations<dim>::energy_component] =
      pressure/(EulerEquations<dim>::gas_gamma - 1.0)
      + 0.5*(u0*u0+v0*v0)*values[EulerEquations<dim>::density_component];
}

//--------------------------------------------------------------------------------------------
// Initial condition for isentropic vortex problem
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void IsentropicVortex<dim>::vector_value (const Point<dim> &p,
                                          Vector<double>   &values) const
{
   const double gamma = EulerEquations<dim>::gas_gamma;
   
   double r2   = (p[0]-x0)*(p[0]-x0) + (p[1]-y0)*(p[1]-y0);
   
   double rho = std::pow(1.0 - a2*exp(1.0-r2), 1.0/(gamma-1.0));
   double vex =  - a1 * (p[1]-y0) * std::exp(0.5*(1.0-r2));
   double vey =  + a1 * (p[0]-x0) * std::exp(0.5*(1.0-r2));
   double pre = std::pow(rho, gamma);
   
   values[0] = rho * vex;
   values[1] = rho * vey;
   values[EulerEquations<dim>::density_component] = rho;
   values[EulerEquations<dim>::energy_component] = pre/(gamma-1.0)
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
   const double gamma = EulerEquations<dim>::gas_gamma;
   
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
   pre = std::pow(rho, gamma);
   
   // Put large pressure in the region { |x| < 0.1 and |y| < 0.1 }
   if(std::fabs(p[0]) < 0.1 && std::fabs(p[1]) < 0.1) pre = 50.0;
   
   values[0] = rho * vex;
   values[1] = rho * vey;
   values[EulerEquations<dim>::density_component] = rho;
   values[EulerEquations<dim>::energy_component] = pre/(gamma-1.0)
                                                 + 0.5 * rho * (vex*vex + vey*vey);
}
//------------------------------------------------------------------------------
double Rayleigh_Taylor_Pressure (double y)
{
   if(y < 0.0)
      return exp(-y/Tl);
   else
      return exp(-y/Tu);
}
//------------------------------------------------------------------------------
double Rayleigh_Taylor_Density (double yc, double y)
{
   if(yc < 0.0)
      return exp(-y/Tl) / Tl;
   else
      return exp(-y/Tu) / Tu;
}
//------------------------------------------------------------------------------
// eta is amplitude of y velocity
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_Rayleigh_Taylor (const double eta)
{
   AssertThrow(dim==2, ExcNotImplemented());
   
   // NOTE: We get multiple sets of same support points since fe is an FESystem
   Quadrature<dim> qsupport (fe.get_unit_support_points());
   FEValues<dim>   fe_values (mapping(), fe, qsupport, update_q_points);
   
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   double rho, pressure, v;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for(; cell != endc; ++cell)
   {
      cell->get_dof_indices(dof_indices);
      fe_values.reinit (cell);
      
      double yc = cell->center()[1];
      
      const std::vector<Point<dim> >& p = fe_values.get_quadrature_points();
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
         const double x = p[i][0];
         const double y = p[i][1];
         unsigned int comp_i = fe.system_to_component_index(i).first;

         switch(comp_i)
         {
            case 0: // x momentum
               old_solution(dof_indices[i]) = 0.0;
               break;
               
            case 1: // y momentum
               rho = Rayleigh_Taylor_Density(yc, y);
               v   = eta * std::sin(4.0*M_PI*x) * std::exp(-50*y*y);
               old_solution(dof_indices[i]) = rho * v;
               break;
               
            case 2: // density
               old_solution(dof_indices[i]) = Rayleigh_Taylor_Density(yc, y);
               break;
               
            case 3: // energy
               pressure = Rayleigh_Taylor_Pressure(y);
               rho = Rayleigh_Taylor_Density(yc, y);
               v   = eta * std::sin(4.0*M_PI*x) * std::exp(-50*y*y);
               old_solution(dof_indices[i]) = pressure/(EulerEquations<dim>::gas_gamma-1) +
                                              0.5 * rho * v * v;
               break;
               
            default:
               AssertThrow(false, ExcMessage("Error in set_initial_condition_Rayleigh_Taylor"));
         }
      }
   }
}

//------------------------------------------------------------------------------
// eta is amplitude of y velocity
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_shocktube ()
{
   AssertThrow(dim==2, ExcNotImplemented());
   
   // NOTE: We get multiple sets of same support points since fe is an FESystem
   Quadrature<dim> qsupport (fe.get_unit_support_points());
   FEValues<dim>   fe_values (mapping(), fe, qsupport, update_q_points);
   
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   double rho, pressure, v;
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for(; cell != endc; ++cell)
   {
      cell->get_dof_indices(dof_indices);
      fe_values.reinit (cell);
      
      double xc = cell->center()[0];
      
      const std::vector<Point<dim> >& p = fe_values.get_quadrature_points();
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
         const double x = p[i][0];
         const double y = p[i][1];
         unsigned int comp_i = fe.system_to_component_index(i).first;
         
         switch(comp_i)
         {
            case 0: // x momentum
            case 1: // y momentum
               old_solution(dof_indices[i]) = 0.0;
               break;
               
            case 2: // density
               if(xc < 0.5)
                  old_solution(dof_indices[i]) = 1.0;
               else
                  old_solution(dof_indices[i]) = 0.125;
               break;
               
            case 3: // energy
               if(xc < 0.5)
                  old_solution(dof_indices[i]) = 2.5;
               else
                  old_solution(dof_indices[i]) = 0.25;
               break;
               
            default:
               AssertThrow(false, ExcMessage("Error in set_initial_condition_Rayleigh_Taylor"));
         }
      }
   }
}

//------------------------------------------------------------------------------
// Sets initial condition based on input file.
// For Qk basis we can just do interpolation.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_Qk ()
{
   if(parameters.ic_function == "rt")
      set_initial_condition_Rayleigh_Taylor ();
//      VectorTools::interpolate(mapping(), dof_handler,
//                               RayleighTaylor<dim>(parameters.gravity), old_solution);
   else if(parameters.ic_function == "rrt")
      VectorTools::interpolate(mapping(), dof_handler,
                               RadialRayleighTaylor<dim>(), old_solution);
   else if(parameters.ic_function == "isohydro")
      VectorTools::interpolate(mapping(), dof_handler,
                               IsothermalHydrostatic<dim>(), old_solution);
   else if(parameters.ic_function == "polyhydro")
      VectorTools::interpolate(mapping(), dof_handler,
                               PolytropicHydrostatic<dim>(1.2), old_solution);
   else if(parameters.ic_function == "isenvort")
      VectorTools::interpolate(mapping(), dof_handler,
                               IsentropicVortex<dim>(5.0, 0.0, 0.0), old_solution);
   else if(parameters.ic_function == "vortsys")
      VectorTools::interpolate(mapping(), dof_handler,
                               VortexSystem<dim>(), old_solution);
   else if(parameters.ic_function == "shocktube")
      set_initial_condition_shocktube ();
   else
      VectorTools::interpolate(mapping(), dof_handler,
                               parameters.initial_conditions, old_solution);
   
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
   std::cout << "Setting initial condition\n";
   set_initial_condition_Qk();
}

template class RayleighTaylor<2>;
template class IsothermalHydrostatic<2>;
template class UnsteadyGravity<2>;
template class IsentropicVortex<2>;
template class VortexSystem<2>;
template class ConservationLaw<2>;
