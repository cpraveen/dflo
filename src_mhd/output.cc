#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

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
void ConservationLaw<dim>::output_results ()
{
   TimerOutput::Scope t(computing_timer, "Output");

   typename MHDEquations<dim>::Postprocessor
   postprocessor (parameters.schlieren_plot);
   
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   
   data_out.add_data_vector (current_solution,
                             MHDEquations<dim>::component_names (),
                             DataOut<dim>::type_dof_data,
                             MHDEquations<dim>::component_interpretation ());
   
   data_out.add_data_vector (current_solution, postprocessor);
   
   Vector<float> subdomain (triangulation.n_active_cells());
   for (unsigned int i=0; i<subdomain.size();++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
   data_out.add_data_vector(subdomain, "subdomain");
   
   data_out.build_patches (mapping(), fe.degree);
   
   static unsigned int output_file_number = 0;
   std::string filename = ("output/solution-" +
                           Utilities::int_to_string(output_file_number,4) +
                           "." +
                           Utilities::int_to_string
                           (triangulation.locally_owned_subdomain(),3));
   std::ofstream outfile ((filename + ".vtu").c_str());
   
   DataOutBase::VtkFlags flags(elapsed_time, output_file_number);
   data_out.set_flags(flags);
   data_out.write_vtu (outfile);
   
   static std::vector< std::vector<std::string> > all_files;
   if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
   {
      std::vector<std::string> filenames;
      for (unsigned int i=0;
           i<Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
         filenames.push_back ("output/solution-" +
                              Utilities::int_to_string (output_file_number, 4) +
                              "." +
                              Utilities::int_to_string (i, 3) +
                              ".vtu");
      all_files.push_back (filenames);
      std::ofstream visit_output ("master_file.visit");
      data_out.write_visit_record(visit_output, all_files);
   }
   
   ++output_file_number;
   
   // Write shock indicator
//   DataOut<dim> shock;
//   shock.attach_dof_handler (dh_cell);
//   shock.add_data_vector (mu_shock, "mu_shock");
//   shock.add_data_vector (shock_indicator, "shock_indicator");
//   shock.build_patches (mapping());
//   if(parameters.output_format == "vtk")     
//   {
//      std::ofstream shock_output ("shock.vtu");
//      shock.write_vtu (shock_output);
//   }
//   else if(parameters.output_format == "tecplot") 
//   {
//      std::ofstream shock_output ("shock.plt");
//      shock.write_tecplot (shock_output);
//   }
   
   // Write cell average solution
//   DataOut<dim> avg;
//   avg.attach_dof_handler (dof_handler0);
//   avg.add_data_vector (cell_average,
//                        EulerEquations<dim>::component_names (),
//                        DataOut<dim>::type_dof_data,
//                        EulerEquations<dim>::component_interpretation ());
//   avg.build_patches (mapping());
//   if(parameters.output_format == "vtk")     
//   {
//      std::ofstream avg_output ("avg.vtu");
//      avg.write_vtu (avg_output);
//   }
//   else if(parameters.output_format == "tecplot") 
//   {
//      std::ofstream avg_output ("avg.plt");
//      avg.write_tecplot (avg_output);
//   }
}

template class ConservationLaw<2>;
