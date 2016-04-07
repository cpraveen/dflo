#ifndef __IC_H__
#define __IC_H__

#include <base/function.h>
#include <lac/vector.h>

#include "equation.h"

//------------------------------------------------------------------------
template <int dim>
class PolytropicHydrostatic: public dealii::Function<dim>
{
public:
   PolytropicHydrostatic (double nu)
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components),
   nu (nu)
   {}
   virtual void vector_value (const dealii::Point<dim>  &p,
                              dealii::Vector<double>  &values) const;
   
private:
   double nu; // polytropic index
   const double p0 = 1.0;
   const double rho0 = 1.0;
   const double alpha = 1.0;
};

//------------------------------------------------------------------------
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
//------------------------------------------------------------------------

template <int dim>
class RadialRayleighTaylor : public dealii::Function<dim>
{
public:
   RadialRayleighTaylor ()
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components)
   {}
   virtual void vector_value (const dealii::Point<dim>  &p,
                              dealii::Vector<double>  &values) const;
   
private:
   const double r0 = 0.6;
   const double eta = 0.02;
   const double k  = 20.0;
   const double drho  = 0.1;
};

//------------------------------------------------------------------------
// Test case from Xing/Shu
//------------------------------------------------------------------------
template <int dim>
class IsothermalHydrostatic : public dealii::Function<dim>
{
public:
   IsothermalHydrostatic (double eta=0.0)
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components),
   eta (eta)
   {}
   virtual void vector_value (const dealii::Point<dim>  &p,
                              dealii::Vector<double>  &values) const;
   
private:
   const double rho0 = 1.21;
   const double p0 = 1.0;
   const double g  = 1.0;
   double eta; // Amplitude of pressure perturbation
};

//------------------------------------------------------------------------
// Test case from Xing/Shu sec. 5.4
//------------------------------------------------------------------------
template <int dim>
class UnsteadyGravity : public dealii::Function<dim>
{
public:
   UnsteadyGravity (double time)
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components),
   time (time)
   {}
   virtual void vector_value (const dealii::Point<dim>  &p,
                              dealii::Vector<double>  &values) const;
   
private:
   double time;
   const double u0 = 1.0;
   const double v0 = 1.0;
   const double p0 = 4.5;
};

//------------------------------------------------------------------------
// Isentropic vortex, just rotates about itself
//------------------------------------------------------------------------
template <int dim>
class IsentropicVortex : public dealii::Function<dim>
{
public:
   IsentropicVortex (double beta, double x0, double y0)
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components),
   beta (beta),
   x0 (x0),
   y0 (y0)
   {
      const double gamma = EulerEquations<dim>::gas_gamma;
      a1 = 0.5*beta/M_PI;
      a2 = (gamma-1.0)*std::pow(a1,2)/2.0;
   }
   virtual void vector_value (const dealii::Point<dim>  &p,
                              dealii::Vector<double>  &values) const;
   
private:
   double beta;
   double a1, a2;
   double x0, y0;
};

//------------------------------------------------------------------------
// Three isentropic vortices
//------------------------------------------------------------------------
template <int dim>
class VortexSystem : public dealii::Function<dim>
{
public:
   VortexSystem ()
   :
   dealii::Function<dim>(EulerEquations<dim>::n_components)
   {
      beta = 5.0;
      Rc = 4.0;
      
      const double gamma = EulerEquations<dim>::gas_gamma;
      a1 = 0.5*beta/M_PI;
      a2 = (gamma-1.0)*std::pow(a1,2)/2.0;
      
      x[0] = 0.0; y[0] = -Rc;
      x[1] = Rc*cos(30.0*M_PI/180.0); y[1] = Rc*sin(30.0*M_PI/180.0);
      x[2] = -x[1]; y[2] = y[1];
   }
   virtual void vector_value (const dealii::Point<dim>  &p,
                              dealii::Vector<double>  &values) const;
   
private:
   double beta, Rc;
   double a1, a2;
   double x[3], y[3];
};

#endif