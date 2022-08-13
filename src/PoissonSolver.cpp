
#include "PoissonSolver.h"
#include "mcmodel.h"

PoissonSolver::PoissonSolver(MeshQuantities * mq) {

  A = mq->init_matrix();

}

Aztec_PoissonSolver::Aztec_PoissonSolver(MeshQuantities * mq) 
  : PoissonSolver(mq) 
{
   Aztec_Problem = new Epetra_LinearProblem;

   Aztec_Problem->SetOperator(A);

}

void Aztec_PoissonSolver::solve_poisson(Epetra_Vector * rhs, Epetra_Vector * sol) {

  Aztec_Problem->SetLHS(sol);

  Aztec_Problem->SetRHS(rhs);

  AztecOO solver(*Aztec_Problem);

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  solver.SetAztecOption(AZ_output, 0);
//  solver.SetAztecOption(AZ_scaling, 2);
  solver.SetAztecOption(AZ_graph_fill, 3);
  solver.SetAztecOption(AZ_overlap, 0);
  solver.SetAztecParam(AZ_ilut_fill, 4.0);
  solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e2);
  int Niters = 320;
  solver.SetAztecOption(AZ_kspace, Niters);
  solver.Iterate(Niters, 1e-14);

}

Aztec_PoissonSolver::~Aztec_PoissonSolver() {
  delete Aztec_Problem;
  delete A;
}

Amesos_PoissonSolver::Amesos_PoissonSolver(MeshQuantities * mq)
  : PoissonSolver(mq) 
{
   Amesos_Problem = new Epetra_LinearProblem;

   Amesos_Problem->SetOperator(A);

   //string SolverType("Amesos_Superludist");
   //string SolverType("Amesos_Lapack");
   
   string SolverType("Amesos_Klu");

   DSolver = Factory.Create(SolverType, *Amesos_Problem);

   assert(DSolver != 0);

   DSolver->SymbolicFactorization();

   DSolver->NumericFactorization();

}

void Amesos_PoissonSolver::solve_poisson(Epetra_Vector * rhs, Epetra_Vector * sol) {

  Amesos_Problem->SetLHS(sol);

  Amesos_Problem->SetRHS(rhs);

  DSolver->Solve();

}

Amesos_PoissonSolver::~Amesos_PoissonSolver() {
  delete A;
  delete DSolver;
  delete Amesos_Problem;
}
