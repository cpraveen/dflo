#include "ic.h"
#include "equation.h"

using namespace dealii;

template <int dim>
void RayleighTaylor<dim>::vector_value (const Point<dim> &p,
                                        Vector<double> values) const
{
   // Density
   if(p[1] < 0.0)
      values[EulerEquations<dim>::density_component] = 1.0;
   else
      values[EulerEquations<dim>::density_component] = 2.0;
   
   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[0] = 0.0;
   
   double vel = 0;
   values[dim-1] = values[EulerEquations<dim>::density_component] * vel;
   
   double pressure = 2.5 - 0.1 * values[EulerEquations<dim>::density_component] * p[1];

   // Energy
   values[EulerEquations<dim>::energy_component] =
      pressure/(EulerEquations<dim>::gas_gamma - 1.0)
      + 0.5 * values[EulerEquations<dim>::density_component] * vel * vel;
}