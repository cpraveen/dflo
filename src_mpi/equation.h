#ifndef __EQUATION_H__
#define __EQUATION_H__

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>

#include <numerics/data_out.h>
#include <numerics/vector_tools.h>
#include <numerics/solution_transfer.h>

#include <lac/parallel_vector.h>

#include <distributed/tria.h>
#include <distributed/grid_refinement.h>
#include <distributed/solution_transfer.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

template <int dim>
struct EulerEquations
{
   // First dim components correspond to momentum
   static const unsigned int n_components             = dim + 2;
   static const unsigned int density_component        = dim;
   static const unsigned int energy_component         = dim+1;
   
   static
   std::vector<std::string>
   component_names ()
   {
      std::vector<std::string> names;
      names.push_back ("XMomentum");
      names.push_back ("YMomentum");
      if(dim==3)
         names.push_back ("ZMomentum");
      names.push_back ("Density");
      names.push_back ("Energy");
      
      return names;
   }
   
   
   static
   std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
   component_interpretation ()
   {
      std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, dealii::DataComponentInterpretation::component_is_part_of_vector);
      data_component_interpretation
      .push_back (dealii::DataComponentInterpretation::component_is_scalar);
      data_component_interpretation
      .push_back (dealii::DataComponentInterpretation::component_is_scalar);
      
      return data_component_interpretation;
   }
   
   // Ratio of specific heats
   static const double gas_gamma;
   
   //---------------------------------------------------------------------------
   // Compute kinetic energy from conserved variables
   //---------------------------------------------------------------------------
   template <typename number, typename InputVector>
   static
   number
   compute_kinetic_energy (const InputVector &W)
   {
      number kinetic_energy = 0;
      for (unsigned int d=0; d<dim; ++d)
         kinetic_energy += *(W.begin()+d) *
                           *(W.begin()+d);
      kinetic_energy *= 0.5/(*(W.begin() + density_component));
      
      return kinetic_energy;
   }
   
   //---------------------------------------------------------------------------
   // Compute pressure from conserved variables
   //---------------------------------------------------------------------------
   template <typename number, typename InputVector>
   static
   number
   compute_pressure (const InputVector &W)
   {
      return ((gas_gamma-1.0) *
              (*(W.begin() + energy_component) -
               compute_kinetic_energy<number>(W)));
   }

   //---------------------------------------------------------------------------
   // Compute maximum eigenvalue
   //---------------------------------------------------------------------------
   template <typename InputVector>
   static
   typename InputVector::value_type
   max_eigenvalue (const InputVector &W)
   {
      typedef typename InputVector::value_type number;

      const number pressure = 
         compute_pressure<number> (W);

      number velocity = 0;
      for (unsigned int d=0; d<dim; ++d)
         velocity += *(W.begin()+d) *
                     *(W.begin()+d);
      velocity = std::sqrt(velocity) / (*(W.begin()+density_component));

      return velocity + std::sqrt(gas_gamma * pressure / (*(W.begin()+density_component)));
   }
   
   //---------------------------------------------------------------------------
   // Compute maximum eigenvalue in normal direction
   //---------------------------------------------------------------------------
   template <typename InputVector>
   static
   typename InputVector::value_type
   max_eigenvalue (const InputVector        &W,
                   const dealii::Point<dim> &normal)
   {
      typedef typename InputVector::value_type number;
      
      const number pressure = compute_pressure<number> (W);
      const number sonic = std::sqrt(gas_gamma * pressure / (*(W.begin()+density_component)));
      
      number velocity = 0;
      for (unsigned int d=0; d<dim; ++d)
         velocity += *(W.begin()+d) * normal[d];
      
      velocity /=  (*(W.begin()+density_component));
      
      return std::fabs(velocity) + sonic;
   }
   
   //---------------------------------------------------------------------------
   // Compute sound speed
   //---------------------------------------------------------------------------
   template <typename InputVector>
   static
   typename InputVector::value_type
   sound_speed (const InputVector &W)
   {
      typedef typename InputVector::value_type number;

      const number pressure = 
         compute_pressure<number> (W);
      return std::sqrt(gas_gamma * pressure / (*(W.begin()+density_component)));
   }
   
   
   //---------------------------------------------------------------------------
   // Compute cartesian components of flux
   //---------------------------------------------------------------------------
   template <typename InputVector, typename number>
   static
   void compute_flux_matrix (const InputVector &W,
                             number            (&flux)[n_components][dim])
   {
      // First compute the pressure that
      // appears in the flux matrix, and
      // then compute the first
      // <code>dim</code> columns of the
      // matrix that correspond to the
      // momentum terms:
      const number pressure = compute_pressure<number> (W);
      
      for (unsigned int d=0; d<dim; ++d)
      {
         for (unsigned int e=0; e<dim; ++e)
            flux[d][e] = W[d] *
                         W[e] /
                         W[density_component];
         
         flux[d][d] += pressure;
      }
      
      // Then the terms for the
      // density (i.e. mass
      // conservation), and,
      // lastly, conservation of
      // energy:
      for (unsigned int d=0; d<dim; ++d)
         flux[density_component][d] = W[d];
      
      for (unsigned int d=0; d<dim; ++d)
         flux[energy_component][d] = W[d] /
                                     W[density_component] *
                                     (W[energy_component] + pressure);
   }
   
   //---------------------------------------------------------------------------
   // Compute flux along normal
   //---------------------------------------------------------------------------
   template <typename InputVector, typename number>
   static
   void normal_flux (const InputVector        &W,
                     const dealii::Point<dim> &normal,
                     number                   (&flux)[n_components])
   {
      const number pressure = compute_pressure<number> (W);
      
      number vdotn = 0.0;
      for (unsigned int d=0; d<dim; ++d)
         vdotn += W[d] * normal[d];
      vdotn /= W[density_component];
      
      flux[density_component] = W[density_component] * vdotn;
      flux[energy_component] = (W[energy_component] + pressure) * vdotn;
      for (unsigned int d=0; d<dim; ++d)
         flux[d] = pressure * normal[d] + W[d] * vdotn;
   }
   
   //---------------------------------------------------------------------------
   // Left and right eigenvector matrices
   // Lx, Rx = along x direction
   // Ly, Ry = along y direction
   // Expressions taken from
   // http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf
   // Note: This is implemented only for 2-D
   //---------------------------------------------------------------------------
   static
   void compute_eigen_matrix (const dealii::Vector<double> &W,
                              double            (&Rx)[n_components][n_components],
                              double            (&Lx)[n_components][n_components],
                              double            (&Ry)[n_components][n_components],
                              double            (&Ly)[n_components][n_components])
   {
      double g1   = gas_gamma - 1.0;
      double rho  = W[density_component];
      double E    = W[energy_component];
      double u    = W[0] / rho;
      double v    = W[1] / rho;
      double q2   = u*u + v*v;
      double p    = g1 * (E - 0.5 * rho * q2);
      double c2   = gas_gamma * p / rho;
      double c    = std::sqrt(c2);
      double beta = 0.5/c2;
      double phi2 = 0.5*g1*q2;
      double h    = c2/g1 + 0.5*q2;
      
      Rx[0][0] = 1;      Rx[0][1] = 0;  Rx[0][2] = 1;     Rx[0][3] = 1;
      Rx[1][0] = u;      Rx[1][1] = 0;  Rx[1][2] = u+c;   Rx[1][3] = u-c;
      Rx[2][0] = v;      Rx[2][1] = -1; Rx[2][2] = v;     Rx[2][3] = v;
      Rx[3][0] = 0.5*q2; Rx[3][1] = -v; Rx[3][2] = h+c*u; Rx[3][3] = h-c*u;
      
      Ry[0][0] = 1;      Ry[0][1] = 0;  Ry[0][2] = 1;     Ry[0][3] = 1;
      Ry[1][0] = u;      Ry[1][1] = 1;  Ry[1][2] = u;     Ry[1][3] = u;
      Ry[2][0] = v;      Ry[2][1] = 0;  Ry[2][2] = v+c;   Ry[2][3] = v-c;
      Ry[3][0] = 0.5*q2; Ry[3][1] = u;  Ry[3][2] = h+c*v; Ry[3][3] = h-c*v;
      
      Lx[0][0] = 1-phi2/c2;       Lx[0][1] = g1*u/c2;       Lx[0][2] = g1*v/c2;    Lx[0][3] = -g1/c2;
      Lx[1][0] = v;               Lx[1][1] = 0;             Lx[1][2] = -1;         Lx[1][3] = 0;
      Lx[2][0] = beta*(phi2-c*u); Lx[2][1] = beta*(c-g1*u); Lx[2][2] = -beta*g1*v; Lx[2][3] = beta*g1;
      Lx[3][0] = beta*(phi2+c*u); Lx[3][1] =-beta*(c+g1*u); Lx[3][2] = -beta*g1*v; Lx[3][3] = beta*g1;
      
      Ly[0][0] = 1-phi2/c2;       Ly[0][1] = g1*u/c2;       Ly[0][2] = g1*v/c2;       Ly[0][3] = -g1/c2;
      Ly[1][0] = -u;              Ly[1][1] = 1;             Ly[1][2] = 0;             Ly[1][3] = 0;
      Ly[2][0] = beta*(phi2-c*v); Ly[2][1] =-beta*g1*u;     Ly[2][2] = beta*(c-g1*v); Ly[2][3] = beta*g1;
      Ly[3][0] = beta*(phi2+c*v); Ly[3][1] =-beta*g1*u;     Ly[3][2] =-beta*(c+g1*v); Ly[3][3] = beta*g1;
      
   }
   
   //---------------------------------------------------------------------------
   // convert from conserved to characteristic variables: W = L*W
   //---------------------------------------------------------------------------
   static
   void transform_to_char (const double           (&L)[n_components][n_components],
                           dealii::Vector<double> &W)
   {
      dealii::Vector<double> V(n_components);
      
      V[0] = W[density_component];
      V[n_components-1] = W[energy_component];
      for(unsigned int d=0; d<dim; ++d)
         V[d+1] = W[d];
      
      W = 0;
      for(unsigned int i=0; i<n_components; ++i)
         for(unsigned int j=0; j<n_components; ++j)
            W[i] += L[i][j] * V[j];
   }
   
   //---------------------------------------------------------------------------
   // convert from characteristic to conserved variables: W = R*W
   //---------------------------------------------------------------------------
   static
   void transform_to_con (const double           (&R)[n_components][n_components],
                          dealii::Vector<double> &W)
   {
      dealii::Vector<double> V(n_components);
      
      V = 0;
      for(unsigned int i=0; i<n_components; ++i)
         for(unsigned int j=0; j<n_components; ++j)
            V[i] += R[i][j] * W[j];

      W[density_component] = V[0];
      W[energy_component] = V[n_components-1];
      for(unsigned int d=0; d<dim; ++d)
         W[d] = V[d+1];
      
   }

   // @sect4{EulerEquations::compute_normal_flux}
   
   // On the boundaries of the
   // domain and across hanging
   // nodes we use a numerical flux
   // function to enforce boundary
   // conditions.  This routine is
   // the basic Lax-Friedrich's flux
   // with a stabilization parameter
   // $\alpha$. It's form has also
   // been given already in the
   // introduction:

   // --------------------------------------------------------------------------
   // Local lax-Friedrichs flux
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void lxf_flux 
   (
    const dealii::Point<dim>         &normal,
    const InputVector                &Wplus,
    const InputVector                &Wminus,
    const dealii::Vector<double>     &Aplus,
    const dealii::Vector<double>     &Aminus,
    typename InputVector::value_type (&normal_flux)[n_components]
   )
   {
      typedef typename InputVector::value_type number;

      // Normal velocity
      number vdotn_plus=0, vdotn_minus=0;
      
      for(unsigned int d=0; d<dim; ++d)
      {
         vdotn_plus  += Wplus[d]  * normal[d];
         vdotn_minus += Wminus[d] * normal[d];
      }
      
      vdotn_plus  /= Wplus [density_component];
      vdotn_minus /= Wminus[density_component];
      
      // pressure
      number p_plus, p_minus;

      p_plus  = compute_pressure<number> (Wplus);
      p_minus = compute_pressure<number> (Wminus);
      
      // Maximum eigenvalue at cell face
      number lambda_plus = max_eigenvalue (Aplus, normal);
      number lambda_minus = max_eigenvalue (Aminus, normal);
      number lambda = std::max(lambda_plus, lambda_minus);
      
      // Momentum flux
      for (unsigned int d=0; d<dim; ++d)
         normal_flux[d] = 0.5 * ( p_plus  * normal[d] + Wplus [d] * vdotn_plus +
                                  p_minus * normal[d] + Wminus[d] * vdotn_minus );

      // Density flux
      normal_flux[density_component] = 0.5 * (Wplus [density_component] * vdotn_plus +
                                              Wminus[density_component] * vdotn_minus);
      
      // Energy flux
      normal_flux[energy_component] = 0.5 * ((Wplus [energy_component] + p_plus)  * vdotn_plus +
                                             (Wminus[energy_component] + p_minus) * vdotn_minus);
      
      // Dissipation flux
      for (unsigned int c=0; c<n_components; ++c)
         normal_flux[c] += 0.5 * lambda * (Wplus[c] - Wminus[c]);
   }
   
   // --------------------------------------------------------------------------
   // Steger-Warming flux
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void steger_warming_flux 
   (
    const dealii::Point<dim>         &normal,
    const InputVector                &Wplus,
    const InputVector                &Wminus,
    typename InputVector::value_type (&normal_flux)[n_components]
   )
   {
      typedef typename InputVector::value_type number;

      number pflux[n_components], mflux[n_components];
      
      // normal velocity and velocity magnitude
      number vdotn_plus=0, vdotn_minus=0, q2_plus=0, q2_minus=0;

      for(unsigned int d=0; d<dim; ++d)
      {
         vdotn_plus  += Wplus[d]  * normal[d];
         vdotn_minus += Wminus[d] * normal[d];
         
         q2_plus  += Wplus[d]  * Wplus[d];
         q2_minus += Wminus[d] * Wminus[d];
      }
      
      vdotn_plus  /= Wplus [density_component];
      vdotn_minus /= Wminus[density_component];
      
      q2_plus  /= Wplus [density_component] * Wplus [density_component];
      q2_minus /= Wminus[density_component] * Wminus[density_component];
      
      // pressure
      number p_plus, p_minus;
      
      p_plus  = compute_pressure<number> (Wplus);
      p_minus = compute_pressure<number> (Wminus);
      
      // sound speed
      number c_plus, c_minus;
      c_plus  = std::sqrt(gas_gamma * p_plus  / Wplus [density_component]);
      c_minus = std::sqrt(gas_gamma * p_minus / Wminus[density_component]);

      // positive flux
      number l1p, l2p, l3p, ap, fp;
      
      l1p = std::max( vdotn_plus,          0.0);
      l2p = std::max( vdotn_plus + c_plus, 0.0);
      l3p = std::max( vdotn_plus - c_plus, 0.0);
      ap  = 2.0 * (gas_gamma - 1.0) * l1p + l2p + l3p;
      fp  = 0.5 * Wplus[density_component] / gas_gamma;
      
      for(unsigned int d=0; d<dim; ++d)
         pflux[d] = ap * Wplus[d]/Wplus[density_component] +
                          c_plus * (l2p - l3p) * normal[d];
      
      pflux[density_component] = ap;
      pflux[energy_component] = 0.5 * ap * q2_plus +
                                c_plus * vdotn_plus * (l2p - l3p) +
                                c_plus * c_plus * (l2p + l3p) / (gas_gamma - 1.0);
      
      // negative flux
      number l1m, l2m, l3m, am, fm;
      
      l1m = std::min( vdotn_minus,           0.0);
      l2m = std::min( vdotn_minus + c_minus, 0.0);
      l3m = std::min( vdotn_minus - c_minus, 0.0);
      am  = 2.0 * (gas_gamma - 1.0) * l1m + l2m + l3m;
      fm  = 0.5 * Wminus[density_component] / gas_gamma;
      
      for(unsigned int d=0; d<dim; ++d)
         mflux[d] = am * Wminus[d]/Wminus[density_component] +
                    c_minus * (l2m - l3m) * normal[d];
      
      mflux[density_component] = am;
      mflux[energy_component] = 0.5 * am * q2_minus +
                                c_minus * vdotn_minus * (l2m - l3m) +
                                c_minus * c_minus * (l2m + l3m) / (gas_gamma - 1.0);
            
      // Total flux
      for (unsigned int c=0; c<n_components; ++c)
         normal_flux[c] = fp * pflux[c] + fm * mflux[c];
   }
   
   // --------------------------------------------------------------------------
   // Roe flux
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void roe_flux
   (
    const dealii::Point<dim>         &normal,
    const InputVector                &W_l,
    const InputVector                &W_r,
    typename InputVector::value_type (&normal_flux)[n_components]
    )
   {
      typedef typename InputVector::value_type number;
      
      number rho_l_sqrt = std::sqrt(W_l[density_component]);
      number rho_r_sqrt = std::sqrt(W_r[density_component]);
      number fact_l = rho_l_sqrt / (rho_l_sqrt + rho_r_sqrt);
      number fact_r = 1.0 - fact_l;
      
      number v_l[dim], v_r[dim], velocity[dim], dv[dim];
      number v2_l = 0, v2_r = 0;
      number v_l_normal = 0, v_r_normal = 0;
      number vel_normal = 0, v2 = 0;
      number v_dot_dv = 0;
      for(unsigned int d=0; d<dim; ++d)
      {
         v_l[d]      = W_l[d] / W_l[density_component];
         v_r[d]      = W_r[d] / W_r[density_component];
         v2_l       += v_l[d] * v_l[d];
         v2_r       += v_r[d] * v_r[d];
         v_l_normal += v_l[d] * normal[d];
         v_r_normal += v_r[d] * normal[d];
         
         velocity[d] = v_l[d] * fact_l + v_r[d] * fact_r;
         vel_normal += velocity[d] * normal[d];
         v2         += velocity[d] * velocity[d];
         dv[d]       = v_r[d] - v_l[d];
         v_dot_dv   += velocity[d] * dv[d];
      }
      
      number p_l = (gas_gamma-1) * (W_l[energy_component] - 0.5 * W_l[density_component] * v2_l);
      number p_r = (gas_gamma-1) * (W_r[energy_component] - 0.5 * W_r[density_component] * v2_r);
      
      number h_l = gas_gamma * p_l / W_l[density_component] / (gas_gamma-1) + 0.5 * v2_l;
      number h_r = gas_gamma * p_r / W_r[density_component] / (gas_gamma-1) + 0.5 * v2_r;
      
      number density = rho_l_sqrt * rho_r_sqrt;
      number h = h_l * fact_l + h_r * fact_r;
      number c = std::sqrt( (gas_gamma-1.0) * (h - 0.5*v2) );
      number drho = W_r[density_component] - W_l[density_component];
      number dp = p_r - p_l;
      number dvn = v_r_normal - v_l_normal;
      
      number a1 = (dp - density * c * dvn) / (2.0*c*c);
      number a2 = drho - dp / (c*c);
      number a3 = (dp + density * c * dvn) / (2.0*c*c);

      number l1 = std::fabs(vel_normal - c);
      number l2 = std::fabs(vel_normal);
      number l3 = std::fabs(vel_normal + c);

      // entropy fix
      number delta = 0.1 * c;
      if(l1 < delta) l1 = 0.5 * (l1*l1/delta + delta);
      if(l3 < delta) l3 = 0.5 * (l3*l3/delta + delta);
      
      number Dflux[n_components];
      Dflux[density_component] = l1 * a1 + l2 * a2 + l3 * a3;
      Dflux[energy_component] = l1 * a1 * (h - c * vel_normal)
                              + l2 * a2 * 0.5 * v2
                              + l2 * density * (v_dot_dv - vel_normal * dvn)
                              + l3 * a3 * (h + c * vel_normal);
      normal_flux[density_component] = 0.5 * (W_l[density_component] * v_l_normal +
                                              W_r[density_component] * v_r_normal
                                              - Dflux[density_component]);
      normal_flux[energy_component] = 0.5 * (W_l[density_component] * h_l * v_l_normal +
                                             W_r[density_component] * h_r * v_r_normal
                                              - Dflux[energy_component]);
      number p_avg = 0.5 * (p_l + p_r);
      for(unsigned int d=0; d<dim; ++d)
      {
         Dflux[d] = (velocity[d] - normal[d] * c) * l1 * a1
                  + velocity[d] * l2 * a2
                  + (dv[d] - normal[d] * dvn) * l2 * density
                  + (velocity[d] + normal[d] * c) * l3 * a3;
         normal_flux[d] = normal[d] * p_avg
                        + 0.5 * (W_l[d] * v_l_normal + W_r[d] * v_r_normal)
                        - 0.5 * Dflux[d];
      }
   }
   
   
   // --------------------------------------------------------------------------
   // HLLC flux
   // Code borrowed from SU2 v2.0.2
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void hllc_flux
   (
    const dealii::Point<dim>         &normal,
    const InputVector                &W_l,
    const InputVector                &W_r,
    typename InputVector::value_type (&normal_flux)[n_components]
    )
   {
      typedef typename InputVector::value_type number;
      
      number rho_l_sqrt = std::sqrt(W_l[density_component]);
      number rho_r_sqrt = std::sqrt(W_r[density_component]);
      number fact_l = rho_l_sqrt / (rho_l_sqrt + rho_r_sqrt);
      number fact_r = 1.0 - fact_l;
      
      number v_l[dim], v_r[dim], velocity[dim];
      number v2_l = 0, v2_r = 0;
      number v_l_normal = 0, v_r_normal = 0;
      number vel_normal = 0, v2 = 0;
      for(unsigned int d=0; d<dim; ++d)
      {
         v_l[d]      = W_l[d] / W_l[density_component];
         v_r[d]      = W_r[d] / W_r[density_component];
         v2_l       += v_l[d] * v_l[d];
         v2_r       += v_r[d] * v_r[d];
         v_l_normal += v_l[d] * normal[d];
         v_r_normal += v_r[d] * normal[d];
         
         velocity[d] = v_l[d] * fact_l + v_r[d] * fact_r;
         vel_normal += velocity[d] * normal[d];
         v2         += velocity[d] * velocity[d];
      }
      
      //pressure
      number p_l = (gas_gamma-1) * (W_l[energy_component] - 0.5 * W_l[density_component] * v2_l);
      number p_r = (gas_gamma-1) * (W_r[energy_component] - 0.5 * W_r[density_component] * v2_r);
      
      // enthalpy
      number h_l = (W_l[energy_component] + p_l) / W_l[density_component];
      number h_r = (W_r[energy_component] + p_r) / W_r[density_component];

      // sound speed
      number c_l = std::sqrt(gas_gamma * p_l / W_l[density_component]);
      number c_r = std::sqrt(gas_gamma * p_r / W_r[density_component]);
      
      // energy per unit mass
      number e_l = W_l[energy_component] / W_l[density_component];
      number e_r = W_r[energy_component] / W_r[density_component];
      
      // roe average
      number h = h_l * fact_l + h_r * fact_r;
      number c = std::sqrt( (gas_gamma-1.0) * (h - 0.5*v2) );
      
      // speed of sound at l and r
      number s_l = std::min(vel_normal-c, v_l_normal-c_l);
      number s_r = std::max(vel_normal+c, v_r_normal+c_r);

      // speed of contact
      number s_m = (p_l - p_r
                    - W_l[density_component] * v_l_normal * (s_l-v_l_normal)
                    + W_r[density_component] * v_r_normal * (s_r-v_r_normal))
      /(W_r[density_component]*(s_r-v_r_normal) - W_l[density_component]*(s_l-v_l_normal));
      
      // Pressure at right and left (Pressure_j=Pressure_i) side of contact surface
      number pStar = W_r[density_component] * (v_r_normal-s_r)*(v_r_normal-s_m) + p_r;

      if (s_m >= 0.0) {
         if (s_l > 0.0)
         {
            normal_flux[density_component] = W_l[density_component]*v_l_normal;
            for (unsigned int d = 0; d < dim; d++)
               normal_flux[d] = W_l[density_component]*v_l[d]*v_l_normal + p_l*normal[d];
            normal_flux[energy_component] = e_l*W_l[density_component]*v_l_normal + p_l*v_l_normal;
         }
         else
         {
            number invSLmSs = 1.0/(s_l-s_m);
            number sLmuL = s_l-v_l_normal;
            number rhoSL = W_l[density_component]*sLmuL*invSLmSs;
            number rhouSL[dim];
            for (unsigned int d = 0; d < dim; d++)
               rhouSL[d] = (W_l[density_component]*v_l[d]*sLmuL+(pStar-p_l)*normal[d])*invSLmSs;
            number eSL = (sLmuL*e_l*W_l[density_component]-p_l*v_l_normal+pStar*s_m)*invSLmSs;
            
            normal_flux[density_component] = rhoSL*s_m;
            for (unsigned int d = 0; d < dim; d++)
               normal_flux[d] = rhouSL[d]*s_m + pStar*normal[d];
            normal_flux[energy_component] = (eSL+pStar)*s_m;
         }
      }
      else
      {
         if (s_r >= 0.0)
         {
            number invSRmSs = 1.0/(s_r-s_m);
            number sRmuR = s_r-v_r_normal;
            number rhoSR = W_r[density_component]*sRmuR*invSRmSs;
            number rhouSR[dim];
            for (unsigned int d = 0; d < dim; d++)
               rhouSR[d] = (W_r[density_component]*v_r[d]*sRmuR+(pStar-p_r)*normal[d])*invSRmSs;
            number eSR = (sRmuR*e_r*W_r[density_component]-p_r*v_r_normal+pStar*s_m)*invSRmSs;
            
            normal_flux[density_component] = rhoSR*s_m;
            for (unsigned int d = 0; d < dim; d++)
               normal_flux[d] = rhouSR[d]*s_m + pStar*normal[d];
            normal_flux[energy_component] = (eSR+pStar)*s_m;
         }
         else
         {
            normal_flux[density_component] = W_r[density_component]*v_r_normal;
            for (unsigned int d = 0; d < dim; d++)
               normal_flux[d] = W_r[density_component]*v_r[d]*v_r_normal + p_r*normal[d];
            normal_flux[energy_component] = e_r*W_r[density_component]*v_r_normal + p_r*v_r_normal;
         }
      }
      
   }
   
   // --------------------------------------------------------------------------
   // Error function
   // --------------------------------------------------------------------------
   template <typename number>
   static
   number ERF(number xarg)
   {
      // constants
      const double a1 =  0.254829592;
      const double a2 = -0.284496736;
      const double a3 =  1.421413741;
      const double a4 = -1.453152027;
      const double a5 =  1.061405429;
      const double p  =  0.3275911;

      // Save the sign of x
      int sign = 1;
      if (xarg < 0)
         sign = -1;
      number x = std::fabs(xarg);

      // A&S formula 7.1.26
      number t = 1.0/(1.0 + p*x);
      number y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

      return sign * y;
   }
   
   // --------------------------------------------------------------------------
   // Kinetic split fluxes
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void kinetic_split_flux
   (
    int                               sign,
    const dealii::Point<dim>         &normal,
    const InputVector                &W,
    typename InputVector::value_type (&normal_flux)[n_components]
   )
   {
      typedef typename InputVector::value_type number;
      
      // normal velocity
      number vdotn=0;

      for(unsigned int d=0; d<dim; ++d)
         vdotn += W[d] * normal[d];
      
      vdotn  /= W[density_component];
      
      // pressure
      number pressure, beta, s, A, B, ufact;
      
      pressure  = compute_pressure<number> (W);
      beta      = 0.5 * W[density_component] / pressure;
      s         = vdotn * std::sqrt(beta);
      A         = 0.5 * (1.0 + sign * ERF(s));
      B         = 0.5 * sign * std::exp(-s*s) / std::sqrt(M_PI * beta);
      ufact     = vdotn * A + B;
      
      for(unsigned int d=0; d<dim; ++d)
         normal_flux[d] = pressure * normal[d] * A + W[d] * ufact;
      
      normal_flux[density_component] = W[density_component] * ufact;
      normal_flux[energy_component]  = (W[energy_component] + pressure) * vdotn * A +
                                       (W[energy_component] + 0.5 * pressure) * B;
      
   }

   // --------------------------------------------------------------------------
   // KFVS flux of Deshpande and Mandal
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void kfvs_flux 
   (
    const dealii::Point<dim>          &normal,
    const InputVector                 &Wplus,
    const InputVector                 &Wminus,
    typename InputVector::value_type  (&normal_flux)[n_components]
   )
   {
      typedef typename InputVector::value_type number;

      number pflux[n_components], mflux[n_components];

      kinetic_split_flux (+1,
                          normal,
                          Wplus,
                          pflux);

      kinetic_split_flux (-1,
                          normal,
                          Wminus,
                          mflux);

      for (unsigned int c=0; c<n_components; ++c)
         normal_flux[c] = pflux[c] + mflux[c];
   }
   
   // --------------------------------------------------------------------------
   // Flux on slip walls. Only pressure flux is present
   // --------------------------------------------------------------------------
   template <typename InputVector>
   static
   void no_penetration_flux 
   (
    const dealii::Point<dim>         &normal,
    const InputVector                &Wminus,
    typename InputVector::value_type (&normal_flux)[n_components]
   )
   {
      typedef typename InputVector::value_type number;

      // pressure
      number pressure = compute_pressure<number> (Wminus);
      
      for (unsigned int c=0; c<n_components; ++c)
         normal_flux[c] = 0.0;
      
      // Only pressure flux is present
      for (unsigned int c=0; c<dim; ++c)
         normal_flux[c] = pressure * normal[c];
   }
   
   //---------------------------------------------------------------------------
   // @sect4{EulerEquations::compute_forcing_vector}
   //---------------------------------------------------------------------------
   
   // In the same way as describing the flux
   // function $\mathbf F(\mathbf w)$, we
   // also need to have a way to describe
   // the right hand side forcing term. As
   // mentioned in the introduction, we
   // consider only gravity here, which
   // leads to the specific form $\mathbf
   // G(\mathbf w) = \left(
   // g_1\rho, g_2\rho, g_3\rho, 0,
   // \rho \mathbf g \cdot \mathbf v
   // \right)^T$, shown here for
   // the 3d case. More specifically, we
   // will consider only $\mathbf
   // g=(0,0,-1)^T$ in 3d, or $\mathbf
   // g=(0,-1)^T$ in 2d. This naturally
   // leads to the following function:
   template <typename InputVector, typename number>
   static
   void compute_forcing_vector (const InputVector &W,
                                number            (&forcing)[n_components])
   {
      const double gravity = -1.0;
      
      for (unsigned int c=0; c<n_components; ++c)
      switch (c)
      {
	      case dim-1:
            forcing[c] = gravity * W[density_component];
            break;

	      case energy_component:
            forcing[c] = gravity * W[dim-1];
            break;

	      default:
            forcing[c] = 0;
      }
   }
   
   
   //---------------------------------------------------------------------------
   // @sect4{Dealing with boundary conditions}
   //---------------------------------------------------------------------------
   
   // Another thing we have to deal with is
   // boundary conditions. To this end, let
   // us first define the kinds of boundary
   // conditions we currently know how to
   // deal with:
   enum BoundaryKind
   {
      inflow_boundary,
      outflow_boundary,
      no_penetration_boundary,
      pressure_boundary,
      farfield_boundary
   };
   
   
   // The next part is to actually decide
   // what to do at each kind of
   // boundary. To this end, remember from
   // the introduction that boundary
   // conditions are specified by choosing a
   // value $\mathbf w^-$ on the outside of
   // a boundary given an inhomogeneity
   // $\mathbf j$ and possibly the
   // solution's value $\mathbf w^+$ on the
   // inside. Both are then passed to the
   // numerical flux $\mathbf
   // H(\mathbf{w}^+, \mathbf{w}^-,
   // \mathbf{n})$ to define boundary
   // contributions to the bilinear form.
   //
   // Boundary conditions can in some cases
   // be specified for each component of the
   // solution vector independently. For
   // example, if component $c$ is marked
   // for inflow, then $w^-_c = j_c$. If it
   // is an outflow, then $w^-_c =
   // w^+_c$. These two simple cases are
   // handled first in the function below.
   //
   // There is a little snag that makes this
   // function unpleasant from a C++
   // language viewpoint: The output vector
   // <code>Wminus</code> will of course be
   // modified, so it shouldn't be a
   // <code>const</code> argument. Yet it is
   // in the implementation below, and needs
   // to be in order to allow the code to
   // compile. The reason is that we call
   // this function at a place where
   // <code>Wminus</code> is of type
   // <code>Table@<2,Sacado::Fad::DFad@<double@>
   // @></code>, this being 2d table with
   // indices representing the quadrature
   // point and the vector component,
   // respectively. We call this function
   // with <code>Wminus[q]</code> as last
   // argument; subscripting a 2d table
   // yields a temporary accessor object
   // representing a 1d vector, just what we
   // want here. The problem is that a
   // temporary accessor object can't be
   // bound to a non-const reference
   // argument of a function, as we would
   // like here, according to the C++ 1998
   // and 2003 standards (something that
   // will be fixed with the next standard
   // in the form of rvalue references).  We
   // get away with making the output
   // argument here a constant because it is
   // the <i>accessor</i> object that's
   // constant, not the table it points to:
   // that one can still be written to. The
   // hack is unpleasant nevertheless
   // because it restricts the kind of data
   // types that may be used as template
   // argument to this function: a regular
   // vector isn't going to do because that
   // one can not be written to when marked
   // <code>const</code>. With no good
   // solution around at the moment, we'll
   // go with the pragmatic, even if not
   // pretty, solution shown here:
   template <typename DataVector>
   static
   void
   compute_Wminus (const BoundaryKind           &boundary_kind,
                   const dealii::Point<dim>     &normal_vector,
                   const DataVector             &Wplus,
                   const dealii::Vector<double> &boundary_values,
                   const DataVector             &Wminus)
   {
      switch (boundary_kind)
      {
	      case inflow_boundary:
	      {
            for (unsigned int c = 0; c < n_components; ++c)
                  Wminus[c] = boundary_values(c);
            break;
	      }
            
	      case outflow_boundary:
	      {
            for (unsigned int c = 0; c < n_components; ++c)
               Wminus[c] = Wplus[c];
            break;
	      }
            
            // Prescribed pressure boundary
            // conditions are a bit more
            // complicated by the fact that
            // even though the pressure is
            // prescribed, we really are
            // setting the energy component
            // here, which will depend on
            // velocity and pressure. So
            // even though this seems like
            // a Dirichlet type boundary
            // condition, we get
            // sensitivities of energy to
            // velocity and density (unless
            // these are also prescribed):
	      case pressure_boundary:
	      {
            const typename DataVector::value_type
            density = Wplus[density_component];
            
            typename DataVector::value_type kinetic_energy = 0;
            for (unsigned int d=0; d<dim; ++d)
                  kinetic_energy += Wplus[d]*Wplus[d];
            kinetic_energy *= 0.5/density;
            
            for (unsigned int c = 0; c < dim; ++c)
               Wminus[c] = Wplus[c];

            Wminus[density_component] = density;
            Wminus[energy_component] = boundary_values(energy_component) / (gas_gamma-1.0) +
                        kinetic_energy;
            
            break;
	      }
            
	      case no_penetration_boundary:
	      {
            // We prescribe the
            // velocity (we are dealing with a
            // particular component here so
            // that the average of the
            // velocities is orthogonal to the
            // surface normal.  This creates
            // sensitivies of across the
            // velocity components.
            typename DataVector::value_type
               vdotn = 0;
            for (unsigned int d = 0; d < dim; d++)
               vdotn += Wplus[d]*normal_vector[d];
            
            for (unsigned int c = 0; c < dim; ++c)
            {
               Wminus[c] = Wplus[c] - 2.0 * vdotn * normal_vector[c];
            }

            Wminus[density_component] = Wplus[density_component];
            Wminus[energy_component]  = Wplus[energy_component];
            break;
	      }
            
	      case farfield_boundary:
	      {
            for (unsigned int c = 0; c < n_components; ++c)
                  Wminus[c] = boundary_values(c);
            break;
	      }
            
	      default:
            Assert (false, dealii::ExcNotImplemented());
      }
   }
   
   
   //---------------------------------------------------------------------------
   // Compute entropy variables V, given conserved variables W
   //---------------------------------------------------------------------------
   template <typename InputVector, typename number>
   static
   void entropy_var (const InputVector &W,
                     number            (&V)[n_components])
   {
      number pressure = compute_pressure<number> (W);
      number T = pressure / W[density_component];

      number u2 = 0;
      for(unsigned int d=0; d<dim; ++d)
      {
         number u = W[d] / W[density_component];
         V[d] = u / T;
         u2 += u * u;
      }

      V[density_component] = log(W[density_component] / std::pow(T, 1.0/(gas_gamma-1.0))) 
                           - 0.5 * u2 / T;
      V[energy_component] = -1.0 / T;
   }
   

   // @sect4{EulerEquations::compute_refinement_indicators}
   
   // In this class, we also want to specify
   // how to refine the mesh. The class
   // <code>ConservationLaw</code> that will
   // use all the information we provide
   // here in the <code>EulerEquation</code>
   // class is pretty agnostic about the
   // particular conservation law it solves:
   // as doesn't even really care how many
   // components a solution vector
   // has. Consequently, it can't know what
   // a reasonable refinement indicator
   // would be. On the other hand, here we
   // do, or at least we can come up with a
   // reasonable choice: we simply look at
   // the gradient of the density, and
   // compute
   // $\eta_K=\log\left(1+|\nabla\rho(x_K)|\right)$,
   // where $x_K$ is the center of cell $K$.
   //
   // There are certainly a number of
   // equally reasonable refinement
   // indicators, but this one does, and it
   // is easy to compute:
   static
   void
   compute_refinement_indicators (const dealii::DoFHandler<dim> 		&dof_handler,
                                  const dealii::Mapping<dim>    		&mapping,
                                  const dealii::parallel::distributed::Vector<double>  	&solution,
                                  dealii::Vector<double>		&refinement_indicators) //dealii::TrilinosWrappers::MPI::Vector
   {
      const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
      std::vector<unsigned int> dofs (dofs_per_cell);
      
      const dealii::QMidpoint<dim>  quadrature_formula;
      const dealii::UpdateFlags update_flags = dealii::update_gradients;
      dealii::FEValues<dim> fe_v (mapping, dof_handler.get_fe(),
                                  quadrature_formula, update_flags);
      
      std::vector<std::vector<dealii::Tensor<1,dim> > >
      dU (1, std::vector<dealii::Tensor<1,dim> >(n_components));
      
      typename dealii::DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      if(cell->is_locally_owned())
      {
         fe_v.reinit(cell);
         fe_v.get_function_gradients (solution, dU);

         refinement_indicators(cell_no)
				  = std::log(1+
                    std::sqrt(dU[0][density_component] *
                              dU[0][density_component]));
      }
      //refinement_indicators.compress(VectorOperation::insert);

   }
   
   
   
   // @sect4{EulerEquations::Postprocessor}
   
   // Finally, we declare a class that
   // implements a postprocessing of data
   // components. The problem this class
   // solves is that the variables in the
   // formulation of the Euler equations we
   // use are in conservative rather than
   // physical form: they are momentum
   // densities $\mathbf m=\rho\mathbf v$,
   // density $\rho$, and energy density
   // $E$. What we would like to also put
   // into our output file are velocities
   // $\mathbf v=\frac{\mathbf m}{\rho}$ and
   // pressure $p=(\gamma-1)(E-\frac{1}{2}
   // \rho |\mathbf v|^2)$.
   //
   // In addition, we would like to add the
   // possibility to generate schlieren
   // plots. Schlieren plots are a way to
   // visualize shocks and other sharp
   // interfaces. The word "schlieren" is a
   // German word that may be translated as
   // "striae" -- it may be simpler to
   // explain it by an example, however:
   // schlieren is what you see when you,
   // for example, pour highly concentrated
   // alcohol, or a transparent saline
   // solution, into water; the two have the
   // same color, but they have different
   // refractive indices and so before they
   // are fully mixed light goes through the
   // mixture along bent rays that lead to
   // brightness variations if you look at
   // it. That's "schlieren". A similar
   // effect happens in compressible flow
   // because the refractive index
   // depends on the pressure (and therefore
   // the density) of the gas.
   //
   // The origin of the word refers to
   // two-dimensional projections of a
   // three-dimensional volume (we see a 2d
   // picture of the 3d fluid). In
   // computational fluid dynamics, we can
   // get an idea of this effect by
   // considering what causes it: density
   // variations. Schlieren plots are
   // therefore produced by plotting
   // $s=|\nabla \rho|^2$; obviously, $s$ is
   // large in shocks and at other highly
   // dynamic places. If so desired by the
   // user (by specifying this in the input
   // file), we would like to generate these
   // schlieren plots in addition to the
   // other derived quantities listed above.
   //
   // The implementation of the algorithms
   // to compute derived quantities from the
   // ones that solve our problem, and to
   // output them into data file, rests on
   // the DataPostprocessor class. It has
   // extensive documentation, and other
   // uses of the class can also be found in
   // step-29. We therefore refrain from
   // extensive comments.
   class Postprocessor : public dealii::DataPostprocessor<dim>
   {
   public:
      Postprocessor (const bool do_schlieren_plot);
      
      virtual
      void
      compute_derived_quantities_vector 
         (const std::vector<dealii::Vector<double> >              &uh,
          const std::vector<std::vector<dealii::Tensor<1,dim> > > &duh,
          const std::vector<std::vector<dealii::Tensor<2,dim> > > &dduh,
          const std::vector<dealii::Point<dim> >                  &normals,
          const std::vector<dealii::Point<dim> >                  &evaluation_points,
          std::vector<dealii::Vector<double> >                    &computed_quantities) const;
      
      virtual std::vector<std::string> get_names () const;
      
      virtual
      std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation () const;
      
      virtual dealii::UpdateFlags get_needed_update_flags () const;
      
      virtual unsigned int n_output_variables() const;
      
   private:
      const bool do_schlieren_plot;
   };
};

#endif
