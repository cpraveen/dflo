#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "equation.h"

// @sect3{Run time parameter handling}

// Our next job is to define a few
// classes that will contain run-time
// parameters (for example solver
// tolerances, number of iterations,
// stabilization parameter, and the
// like). One could do this in the
// main class, but we separate it
// from that one to make the program
// more modular and easier to read:
// Everything that has to do with
// run-time parameters will be in the
// following namespace, whereas the
// program logic is in the main
// class.
//
// We will split the run-time
// parameters into a few separate
// structures, which we will all put
// into a namespace
// <code>Parameters</code>. Of these
// classes, there are a few that
// group the parameters for
// individual groups, such as for
// solvers, mesh refinement, or
// output. Each of these classes have
// functions
// <code>declare_parameters()</code>
// and
// <code>parse_parameters()</code>
// that declare parameter subsections
// and entries in a ParameterHandler
// object, and retrieve actual
// parameter values from such an
// object, respectively. These
// classes declare all their
// parameters in subsections of the
// ParameterHandler.
//
// The final class of the following
// namespace combines all the
// previous classes by deriving from
// them and taking care of a few more
// entries at the top level of the
// input file, as well as a few odd
// other entries in subsections that
// are too short to warrent a
// structure by themselves.
//
// It is worth pointing out one thing here:
// None of the classes below have a
// constructor that would initialize the
// various member variables. This isn't a
// problem, however, since we will read all
// variables declared in these classes from
// the input file (or indirectly: a
// ParameterHandler object will read it from
// there, and we will get the values from
// this object), and they will be initialized
// this way. In case a certain variable is
// not specified at all in the input file,
// this isn't a problem either: The
// ParameterHandler class will in this case
// simply take the default value that was
// specified when declaring an entry in the
// <code>declare_parameters()</code>
// functions of the classes below.
namespace Parameters
{
   
   // @sect4{Parameters::Solver}
   //
   // The first of these classes deals
   // with parameters for the linear
   // inner solver. It offers
   // parameters that indicate which
   // solver to use (GMRES as a solver
   // for general non-symmetric
   // indefinite systems, or a sparse
   // direct solver), the amount of
   // output to be produced, as well
   // as various parameters that tweak
   // the thresholded incomplete LU
   // decomposition (ILUT) that we use
   // as a preconditioner for GMRES.
   //
   // In particular, the ILUT takes
   // the following parameters:
   // - ilut_fill: the number of extra
   //   entries to add when forming the ILU
   //   decomposition
   // - ilut_atol, ilut_rtol: When
   //   forming the preconditioner, for
   //   certain problems bad conditioning
   //   (or just bad luck) can cause the
   //   preconditioner to be very poorly
   //   conditioned.  Hence it can help to
   //   add diagonal perturbations to the
   //   original matrix and form the
   //   preconditioner for this slightly
   //   better matrix.  ATOL is an absolute
   //   perturbation that is added to the
   //   diagonal before forming the prec,
   //   and RTOL is a scaling factor $rtol
   //   \geq 1$.
   // - ilut_drop: The ILUT will
   //   drop any values that
   //   have magnitude less than this value.
   //   This is a way to manage the amount
   //   of memory used by this
   //   preconditioner.
   //
   // The meaning of each parameter is
   // also briefly described in the
   // third argument of the
   // ParameterHandler::declare_entry
   // call in
   // <code>declare_parameters()</code>.
   struct Solver
   {
      enum SolverType { gmres, direct, umfpack, rk3, mood };
      SolverType solver;
      bool implicit;
      
      enum  OutputType { quiet, verbose };
      OutputType output;
      
      double linear_residual;
      int max_iterations;
      
      double ilut_fill;
      double ilut_atol;
      double ilut_rtol;
      double ilut_drop;
      
      static void declare_parameters (dealii::ParameterHandler &prm);
      void parse_parameters (dealii::ParameterHandler &prm);
   };
   
   
   
   // @sect4{Parameters::Refinement}
   //
   // Similarly, here are a few parameters
   // that determine how the mesh is to be
   // refined (and if it is to be refined at
   // all). For what exactly the shock
   // parameters do, see the mesh refinement
   // functions further down.
   struct Refinement
   {
      bool do_refine;
      double shock_val;
      double shock_levels;
      double refine_time_step;
      int refine_iter_step;
      
      static void declare_parameters (dealii::ParameterHandler &prm);
      void parse_parameters (dealii::ParameterHandler &prm);
   };
   
   
   
   // @sect4{Parameters::Flux}
   //
   // Next a section on flux modifications to
   // make it more stable. In particular, two
   // options are offered to stabilize the
   // Lax-Friedrichs flux: either choose
   // $\mathbf{H}(\mathbf{a},\mathbf{b},\mathbf{n})
   // =
   // \frac{1}{2}(\mathbf{F}(\mathbf{a})\cdot
   // \mathbf{n} + \mathbf{F}(\mathbf{b})\cdot
   // \mathbf{n} + \alpha (\mathbf{a} -
   // \mathbf{b}))$ where $\alpha$ is either a
   // fixed number specified in the input
   // file, or where $\alpha$ is a mesh
   // dependent value. In the latter case, it
   // is chosen as $\frac{h}{2\delta T}$ with
   // $h$ the diameter of the face to which
   // the flux is applied, and $\delta T$
   // the current time step.
   struct Flux
   {
      enum FluxType {lxf, sw, kfvs, roe, hllc, kep};
      FluxType flux_type;

      enum StabilizationKind { constant, mesh_dependent };
      StabilizationKind stabilization_kind;
      
      double stabilization_value;
      
      static void declare_parameters (dealii::ParameterHandler &prm);
      void parse_parameters (dealii::ParameterHandler &prm);
   };
   
   struct Limiter
   {
      enum LimiterType { none, TVB, minmax };
      enum ShockIndType { limiter, density, energy, u2 };
      
      LimiterType  limiter_type;
      ShockIndType shock_indicator_type;
      bool         char_lim;
      bool         pos_lim;
      double       M;
      double       beta;
      bool         conserve_angular_momentum;
      
      static void declare_parameters (dealii::ParameterHandler &prm);
      void parse_parameters (dealii::ParameterHandler &prm);
   };
   
   
   // @sect4{Parameters::Output}
   //
   // Then a section on output parameters. We
   // offer to produce Schlieren plots (the
   // squared gradient of the density, a tool
   // to visualize shock fronts), and a time
   // interval between graphical output in
   // case we don't want an output file every
   // time step.
   struct Output
   {
      bool schlieren_plot;
      double output_time_step;
      int output_iter_step;
      std::string output_format;
      unsigned int ang_mom_step;
      
      static void declare_parameters (dealii::ParameterHandler &prm);
      void parse_parameters (dealii::ParameterHandler &prm);
   };
   
   
   
   // @sect4{Parameters::AllParameters}
   //
   // Finally the class that brings it all
   // together. It declares a number of
   // parameters itself, mostly ones at the
   // top level of the parameter file as well
   // as several in section too small to
   // warrant their own classes. It also
   // contains everything that is actually
   // space dimension dependent, like initial
   // or boundary conditions.
   //
   // Since this class is derived from all the
   // ones above, the
   // <code>declare_parameters()</code> and
   // <code>parse_parameters()</code>
   // functions call the respective functions
   // of the base classes as well.
   //
   // Note that this class also handles the
   // declaration of initial and boundary
   // conditions specified in the input
   // file. To this end, in both cases,
   // there are entries like "w_0 value"
   // which represent an expression in terms
   // of $x,y,z$ that describe the initial
   // or boundary condition as a formula
   // that will later be parsed by the
   // FunctionParser class. Similar
   // expressions exist for "w_1", "w_2",
   // etc, denoting the <code>dim+2</code>
   // conserved variables of the Euler
   // system. Similarly, we allow up to
   // <code>max_n_boundaries</code> boundary
   // indicators to be used in the input
   // file, and each of these boundary
   // indicators can be associated with an
   // inflow, outflow, or pressure boundary
   // condition, with inhomogenous boundary
   // conditions being specified for each
   // component and each boundary indicator
   // separately.
   //
   // The data structure used to store the
   // boundary indicators is a bit
   // complicated. It is an array of
   // <code>max_n_boundaries</code> elements
   // indicating the range of boundary
   // indicators that will be accepted. For
   // each entry in this array, we store a
   // pair of data in the
   // <code>BoundaryCondition</code>
   // structure: first, an array of size
   // <code>n_components</code> that for
   // each component of the solution vector
   // indicates whether it is an inflow,
   // outflow, or other kind of boundary,
   // and second a FunctionParser object
   // that describes all components of the
   // solution vector for this boundary id
   // at once.
   //
   // The <code>BoundaryCondition</code>
   // structure requires a constructor since
   // we need to tell the function parser
   // object at construction time how many
   // vector components it is to
   // describe. This initialization can
   // therefore not wait till we actually
   // set the formulas the FunctionParser
   // object represents later in
   // <code>AllParameters::parse_parameters()</code>
   //
   // For the same reason of having to tell
   // Function objects their vector size at
   // construction time, we have to have a
   // constructor of the
   // <code>AllParameters</code> class that
   // at least initializes the other
   // FunctionParser object, i.e. the one
   // describing initial conditions.
   template <int dim>
   struct AllParameters : public Solver,
   public Refinement,
   public Flux,
   public Limiter,
   public Output
   {
      static const unsigned int max_n_boundaries = 10;
      
      struct BoundaryConditions
      {
         typename EulerEquations<dim>::BoundaryKind kind;
         
         dealii::FunctionParser<dim> values;
         
         BoundaryConditions ();
      };
      
      
      AllParameters ();
      
      double diffusion_power;
      double diffusion_coef;
      
      double gravity;
      dealii::FunctionParser<dim> external_force;
      
      unsigned int degree;
      enum BasisType { Qk, Pk };
      BasisType basis;
      enum MappingType { q1, q2, cartesian };
      MappingType mapping_type;
      
      double cfl;
      std::string time_step_type;
      double time_step, final_time;
      double theta;
      bool is_stationary;
      unsigned int max_nonlin_iter;
      
      std::string mesh_type;
      std::string mesh_filename;
      
      std::string ic_function;
      dealii::FunctionParser<dim> initial_conditions;
      BoundaryConditions  boundary_conditions[max_n_boundaries];
      
      static void declare_parameters (dealii::ParameterHandler &prm);
      void parse_parameters (dealii::ParameterHandler &prm);
   };
   
   
}

#endif
