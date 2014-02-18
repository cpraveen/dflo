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

#include <Sacado.hpp>

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
      velocity = std::sqrt(velocity) / (*W.begin()+density_component);

      return velocity + std::sqrt(gas_gamma * pressure / (*W.begin()+density_component));
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
      return std::sqrt(gas_gamma * pressure / (*W.begin()+density_component));
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
      
      // Maximum eigenvalue at cell face, based on an average
      number vdotn_avg, p_avg, d_avg, lambda;
      vdotn_avg = 0.5 * (vdotn_plus + vdotn_minus);
      p_avg     = 0.5 * (p_plus + p_minus);
      d_avg     = 0.5 * (Wplus[density_component] + Wminus[density_component]);

      lambda = std::fabs(vdotn_avg) + 
               std::sqrt( gas_gamma * p_avg / d_avg );
            
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
      std::cout << "out \n";
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
      //const double gravity = -1.0;
      const double gravity = 0.0;

      
      for (unsigned int c=0; c<n_components; ++c)
      switch (c)
      {
	      case dim-1:
            forcing[c] = gravity * W[density_component];
            break;

	      case energy_component:
            forcing[c] = gravity *
                         W[density_component] *
                         W[dim-1];
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
               Wminus[c] = Wplus[c] - 2.0 * vdotn * normal_vector[c];

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
   compute_refinement_indicators (const dealii::DoFHandler<dim> &dof_handler,
                                  const dealii::Mapping<dim>    &mapping,
                                  const dealii::Vector<double>  &solution,
                                  dealii::Vector<double>        &refinement_indicators)
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
      {
         fe_v.reinit(cell);
         fe_v.get_function_grads (solution, dU);
         
         refinement_indicators(cell_no)
	      = std::log(1+
                    std::sqrt(dU[0][density_component] *
                              dU[0][density_component]));
      }
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
