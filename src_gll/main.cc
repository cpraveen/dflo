#include "claw.h"
#include <base/timer.h>

using namespace dealii;

// @sect3{main()}
// The following ``main'' function is
// similar to previous examples and
// need not to be commented on. Note
// that the program aborts if no input
// file name is given on the command
// line.
int main (int argc, char *argv[])
{
   deallog.depth_console(0);
   if (argc != 2 && argc != 3)
   {
      std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
      std::cout << "      " << argv[0] << " input_file  num_threads" << std::endl;
      std::exit(1);
   }
   unsigned int n_threads = 1;
   if(argc==3) n_threads = atoi(argv[2]);
   
   try
   {
      Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv, n_threads);
      ParameterHandler prm;
      Parameters::AllParameters<2>::declare_parameters (prm);
      bool status = prm.read_input (argv[1], true);
      AssertThrow( status, ExcFileNotOpen(argv[1]) );
      prm.print_parameters(std::cout, ParameterHandler::Text);
      unsigned int degree  = prm.get_integer("degree"); // Degree of FEM
      std::string basis = prm.get("basis");
      
      Timer timer;
      timer.start ();
      
      const FE_DGQArbitraryNodes<2> fe_scalar(QGaussLobatto<1>(degree+1));
      ConservationLaw<2> cons (argv[1], degree, fe_scalar);
      cons.run ();
      //cons.compute_errors ();
      timer.stop ();

      std::cout << std::endl;
      std::cout << "Elapsed CPU time : " << timer()/60 << " min.\n";
      std::cout << "Elapsed wall time: " << timer.wall_time()/60 << " min.\n";
   }
   catch (std::exception &exc)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
   };
   
   return 0;
}
