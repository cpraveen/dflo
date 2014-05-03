#ifndef __IC_H__
#define __IC_H__

#include <base/function.h>
#include <lac/vector.h>

#include "equation.h"

template <int dim>
class RayleighTaylor : public dealii::Function<dim>
{
public:
   RayleighTaylor (double gravity)
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components),
   gravity (gravity)
   {}
   virtual void vector_value (const dealii::Point<dim>  &p,
                                dealii::Vector<double>  &values) const;
   
private:
   double gravity;
   const double Lx = 0.5;  // size of domain in x
   const double Ly = 1.5;  // size of domain in y
   const double A  = 0.01; // y velocity perturbation amplitude
   const double P0 = 2.5;  // pressure at y=0
};

#endif