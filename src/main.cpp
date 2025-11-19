#include "mcmodel.h"

int main(int argc, char *argv[])
{
  
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
  mq->run();
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

return 0;
}
