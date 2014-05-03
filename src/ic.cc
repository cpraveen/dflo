#include "ic.h"
#include "equation.h"

using namespace dealii;

// Initial condition for Rayleigh-Taylor problem
// This is setup for 2-d case only
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

template class RayleighTaylor<2>;