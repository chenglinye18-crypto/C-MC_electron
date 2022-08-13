
#include "utils.h"

class MeshQuantities;

class PoissonSolver {
  public:
    virtual void solve_poisson(Epetra_Vector *, Epetra_Vector *) = 0;

    PoissonSolver(MeshQuantities *);

    virtual ~PoissonSolver(){
      delete A;
    }

  protected: 
    Epetra_CrsMatrix * A;
};

class Aztec_PoissonSolver : public PoissonSolver {
  public:
    Aztec_PoissonSolver(MeshQuantities * mq);

    virtual void solve_poisson(Epetra_Vector *, Epetra_Vector *);

    ~Aztec_PoissonSolver();

  private:

    Epetra_LinearProblem * Aztec_Problem;

};

class Amesos_PoissonSolver : public PoissonSolver {

  public:

    Amesos_PoissonSolver(MeshQuantities * mq);

    virtual void solve_poisson(Epetra_Vector *, Epetra_Vector *);

    ~Amesos_PoissonSolver();

  private:
    Amesos_BaseSolver * DSolver;

    Amesos Factory;

    Epetra_LinearProblem * Amesos_Problem;

};
