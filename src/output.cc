#include <base/utilities.h>

#include <dofs/dof_handler.h>

#include <numerics/data_out.h>

#include <iostream>
#include <fstream>

#include "equation.h"
#include "claw.h"

using namespace dealii;

// @sect4{ConservationLaw::output_results}

// This function now is rather
// straightforward. All the magic, including
// transforming data from conservative
// variables to physical ones has been
// abstracted and moved into the
// EulerEquations class so that it can be
// replaced in case we want to solve some
// other hyperbolic conservation law.
//
// Note that the number of the output file is
// determined by keeping a counter in the
// form of a static variable that is set to
// zero the first time we come to this
// function and is incremented by one at the
// end of each invokation.
template <int dim>
void ConservationLaw<dim>::output_results () const
{
   typename EulerEquations<dim>::Postprocessor
   postprocessor (parameters.schlieren_plot);
   
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   
   data_out.add_data_vector (current_solution,
                             EulerEquations<dim>::component_names (),
                             DataOut<dim>::type_dof_data,
                             EulerEquations<dim>::component_interpretation ());
   
   data_out.add_data_vector (current_solution, postprocessor);
   
   data_out.build_patches (fe.degree);
   
   static unsigned int output_file_number = 0;
   std::string filename = "solution-" + Utilities::int_to_string (output_file_number, 3);
   
   if(parameters.output_format == "vtk")     
      filename += ".vtk";
   else if(parameters.output_format == "tecplot") 
      filename += ".plt";
   
   std::ofstream output (filename.c_str());
   
   if(parameters.output_format == "vtk")     
      data_out.write_vtk (output);
   else if(parameters.output_format == "tecplot") 
      data_out.write_tecplot (output);
   
   ++output_file_number;
   
   // Write shock indicator
   DataOut<dim> shock;
   shock.attach_dof_handler (dh_cell);
   shock.add_data_vector (mu_shock, "shock");
   shock.build_patches ();
   std::ofstream shock_output ("shock.plt");
   shock.write_tecplot (shock_output);
   
   // Write cell average solution
   DataOut<dim> avg;
   avg.attach_dof_handler (dof_handler0);
   avg.add_data_vector (cell_average,
                        EulerEquations<dim>::component_names (),
                        DataOut<dim>::type_dof_data,
                        EulerEquations<dim>::component_interpretation ());
   avg.build_patches ();
   std::ofstream avg_output ("avg.plt");
   avg.write_tecplot (avg_output);
}

template class ConservationLaw<2>;
