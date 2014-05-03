#ifndef __IC_H__
#define __IC_H__

#include <base/function.h>
#include <lac/vector.h>

template <int dim>
class RayleighTaylor : public dealii::Function<dim>
{
public:
   RayleighTaylor () : dealii::Function<dim>() {}
   virtual void vector_value (const dealii::Point<dim>   &p,
                                dealii::Vector<double> values) const;
};

#endif