#include "mcmodel.h"

int main(int argc, char *argv[])
{

  cout << "----------------------------------------------------" << endl;
  cout << "-3D Monte Carlo Simulator for Semiconductor Devices-" << endl;
  cout << "-    Developed by Chenglin Ye, Peking University   -" << endl;
  cout << "-    Email:chenglinye18@gmail.com                  -" << endl;
  cout << "-    Version: 1.0.0, 2025-11-01                    -" << endl;
  cout << "----------------------------------------------------" << endl;
  
#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  Epetra_SerialComm Comm;

#endif

  /* initialize the random generator */
  mysrand(time(NULL));

  int MyPID = Comm.MyPID();
  bool verbose = (MyPID==0);

  

  // Get the number of local equations from the command line
  if (argc!=2)
   {
     if (verbose) 
       cout << "input file missing, Usage: " << argv[0] << " input_file" << endl;
    std::exit(1);
   }

  /*define a object*/
  MeshQuantities * mq = new MeshQuantities;

  /*initialize first*/
  mq->initialize(argv[1]);

  /*begin to run*/
  if (verbose)
    cout << "Calling run()" << endl;
  mq->run();
  if (verbose)
    cout << "run() finished" << endl;
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

return 0;
}

