#ifndef __UTILS_HEADER__
#define __UTILS_HEADER__

#include <sstream>
#include <vector>
#include <list>

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"


#ifndef __cplusplus
#define __cplusplus
#endif

#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"

#include "Trilinos_Util.h"
#include "Trilinos_Util_CommandLineParser.h"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_ConfigDefs.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"

#include "defs.h"
#include "Constant_GlobalVariable_MathFunction.h"

using namespace std;

/**
 * c_numx: Number of Cell in a row in the X directions
 * c_numy: Number of Cell in a row in the Y directions
 * c_numz: Number of Cell in a row in the Z directions 
 * 
 * p_numx: Number of Points in a row in the X directions
 * p_numy: Number of Points in a row in the Y directions
 * p_numz: Number of Points in a row in the Z directions
 */
extern int c_numx, c_numy, c_numz, p_numx, p_numy, p_numz, c_numxz, p_numxz;

extern int c_jbegin, c_jend, c_ibegin, c_isize, c_iend, c_kbegin, c_kend, c_ksize;
extern int p_jbegin, p_jend, p_isize, p_ibegin, p_iend, p_kbegin, p_kend, p_ksize;

extern int c_jbegin_ghost, c_jend_ghost, c_num_local_y_ghost, c_num_local_ghost, p_num_local_ghost, p_jbegin_ghost, p_jend_ghost;

extern int p_jbegin_nonoverlap, p_num_local_nonoverlap, p_num_local_nonoverlap_y, p_jend_nonoverlap;

extern int p_qc_jbegin, p_qc_jend, p_qc_numy;

extern int c_num_local_y, p_num_local_y, c_num_local, p_num_local;

extern int p_tox, p_gbegin, p_gend, p_box;

extern vector<double>  dx, dy, dz;
  
extern vector<double> lx, ly, lz;

extern Epetra_Map * p_map, * c_map, *c_map_y, *c_map_ghost, * p_map_ghost, * c_map_ghost_4, * p_qc_map;

extern Epetra_Map * p_map_nonoverlap;

extern double eps[MNMaterialType];

extern double Vg, Vs, Vd, Eg, Eb;

extern double fermi[2];

extern double tsi, fermi_order;

extern  int mpi_size, mpi_rank;

extern vector<double> int_val;

extern double qc_xratio , qc_xtheta, qc_yratio , qc_ytheta, qc_zratio, qc_ztheta;

extern "C" {
  double dspevx_(char *jobz, char *range, char *uplo, int *n,
                 double *ap, double *vl, double *vu, int *il, int*
                 iu, double *abstol, int *m, double *w, double *z__,
                 int *ldz, double *work, int *iwork, int *ifail,
                 int *info);
  void fermid_(double *, double *, double *, double *, int *);
}

class user_cmd {
public:
  int type;
  int range[6];
  int int_param[10];
  double dbl_param[10];
};

enum {
  WRONG_CELL,
  TOO_MANY_LOOPS,
  EIG_CRASH,
  FAIL_CONVERGE,
  WRONG_CONTACT,
  WRONG_RANDOM,
  GHOST_PARTICLE,
};

class Particle{
public:
  double x,y,z;
  double kx,ky,kz;
  double charge, energy;
  double left_time;
  int i, j,k;
  int itet, isym, par_type;
  int flag;
  int seed;
  int par_id;
};


extern Epetra_CrsMatrix * init_matrix(bool float_boundary);
extern void print_p_data(Epetra_Vector *, string, double scale = 1);
extern void print_c_data(Epetra_Vector *, string, bool, double scale = 1);
extern  void print_c_data(Epetra_IntVector *, string, bool, double scale = 1);


extern double fermi_integral(double x, double order, bool aprox);
extern void output_fermi_integral(double , double);
extern double anti_dummy(double x, double order);
extern double integral(double zeta_s, double zeta_d, double upper_lim);

extern void init_int_val();
extern double find_int_val(double x);

extern bool between(int i, int j, int k);

extern  int sign(double x);

extern string getFileName(const string &, int);

extern void err_message(int, string func_name);

extern Particle make_particle(double * , int * );

extern void cp_to_buf(double * , int * , Particle ) ;

extern double ran1(int *);
extern double Random();
extern void mysrand(int);
extern int ran_seed[1];
extern double fermi_fun(double x);
extern double random_fermi_fun(double Ef, double Emax);
extern double erf(double);

inline double Random(int & seed){
  int Big = 2147483647;
  long long a = 16807;
  long long tmp = a * seed;
  seed = tmp % Big;
  return seed * 1.0 / Big;
}

#endif
