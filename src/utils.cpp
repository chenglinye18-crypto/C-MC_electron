#include "utils.h"

int c_numx, c_numy, c_numz, p_numx, p_numy, p_numz, c_numxz, p_numxz;

int c_jbegin, c_jend, c_ibegin, c_isize, c_iend, p_jbegin, p_jend, p_isize, p_ibegin, p_iend;
int c_kbegin, c_kend, c_ksize, p_kbegin, p_kend, p_ksize;

int c_jbegin_ghost, c_jend_ghost, c_num_local_y_ghost, c_num_local_ghost, p_num_local_ghost, p_jbegin_ghost, p_jend_ghost;

int p_jbegin_nonoverlap, p_num_local_nonoverlap, p_num_local_nonoverlap_y, p_jend_nonoverlap, p_qc_jbegin, p_qc_jend, p_qc_numy;

int c_num_local_y, p_num_local_y, c_num_local, p_num_local;

int p_tox, p_gbegin, p_gend, p_box;

vector<double>  dx, dy, dz;
  
vector<double> lx, ly, lz;

Epetra_Map * p_map, * c_map, *c_map_y, *c_map_ghost, * p_map_ghost, * c_map_ghost_4, *p_qc_map;

Epetra_Map * p_map_nonoverlap;

double eps[MNMaterialType];

double Eg, Vg, Vs, Vd, Eb;

double fermi[2];

double tsi, fermi_order;

int size, rank;

double qc_xratio , qc_xtheta, qc_yratio , qc_ytheta, qc_zratio, qc_ztheta;

vector<double> int_val;


void err_message(int errid, string err_info)
{
  cout << STARS << "ERROR in function: " << err_info << STARS << endl;
  switch (errid) {
    case WRONG_CELL:
      cout << " wrong cell\n";
      break;
    case TOO_MANY_LOOPS:
      cout << " too many loops\n";
      break;
    case EIG_CRASH:
      cout << " eigvalue solver crashed\n";
      break;
    case FAIL_CONVERGE:
      cout << " fail to converge\n";
      break;
    case GHOST_PARTICLE:
      cout << err_info << endl;
      break;
    case WRONG_RANDOM:
      cout << "random should bewteen [0,1] " << err_info << endl;
      break;
  }
}

 int sign(double x)
{
    return int(fabs(x)/x);
}

bool between(int i,int j, int k)
{
  if ((i >= j) && (i <= k)) return true; else return false;
}

string getFileName(const string & prefix, const int num){
  string st_num;
  stringstream ss;
  ss << num;
  ss >> st_num;
  return prefix + st_num;
}

/* for debug use, output the fermi dirac integral of order 0 and 1/2, using two different
 * ways */
void output_fermi_integral(double b, double e) {
  int N = 1000; 
  double h = (e - b) / N;
  double x;
  int i;
  ofstream ofile;
  string filename("fermi_integral");
  filename = "./data/" + filename;
  ofile.open(filename.c_str());

  ofile << N << endl;
  for (i = 0;i < N ;i ++) {
    x = b + i * h; 
    ofile << x << ' ' << fermi_integral(x, 0, true) << ' ' << fermi_integral(x, 0.5, true) << ' ' << fermi_integral(x, 0, false) << ' ' << fermi_integral(x, 0.5, false) << endl;
  }
  ofile.close();
}

double integral(double zeta_s, double zeta_d, double upper_lim) {
  double E_tail = 12, dx, y1, y2, y3;
  double tem_lim, tem_lim_l, tem_lim_h;
  int nx;
  vector<double> x, dummy1;
  int i;
  if (upper_lim > 0) {
    tem_lim = pow(upper_lim , 0.5);
    dx = 1e-2;
    nx = int (tem_lim / dx) + 2;
    dx = tem_lim / (nx - 1);
    x.resize(nx);
    dummy1.resize(nx);
    for (i = 0;i < nx; i ++)
      x[i] = dx * i;
    for (i = 0;i < nx; i ++)
      dummy1[i] = fermi_integral(zeta_s - x[i] * x[i], -0.5, 1); 
    y1 = 0;
    for (i = 0;i < nx; i ++)
      y1 += dummy1[i];

    y1 -= (dummy1[0] + dummy1[nx - 1]) / 2;
    y1 *= dx;
  }

  tem_lim_l = pow(upper_lim, 0.5) ; 
  tem_lim_h = pow(max(upper_lim + E_tail , zeta_d + E_tail), 0.5);
  dx = 1e-2;
  nx = int((tem_lim_h - tem_lim_l) / dx) + 2;
  dx = (tem_lim_h - tem_lim_l) / (nx - 1);
  x.resize(nx);
  dummy1.resize(nx);
  for (i = 0;i < nx; i ++)
    x[i] = tem_lim_l + i * dx;
  dx = x[1] - x[0];
  for (i = 0;i < nx; i ++)
    dummy1[i] = fermi_integral(zeta_d - x[i] * x[i], -0.5, 1); 
  y2 = 0;
  for (i = 0;i < nx; i ++)
    y2 += dummy1[i];

  y2 -= (dummy1[0] + dummy1[nx - 1]) / 2;
  y2 *= dx;

  tem_lim_l = 0;
  tem_lim_h = pow(max(upper_lim + E_tail , zeta_s + E_tail), 0.5);
  dx = 1e-2;
  nx = int((tem_lim_h - tem_lim_l) / dx) + 2;
  dx = (tem_lim_h - tem_lim_l) / (nx - 1);
  x.resize(nx);
  dummy1.resize(nx);
  for (i = 0;i < nx; i ++)
    x[i] = tem_lim_l + i * dx;
  dx = x[1] - x[0];
  for (i = 0;i < nx; i ++)
    dummy1[i] = fermi_integral(zeta_s - x[i] * x[i], -0.5, 1); 
  y3 = 0;
  for (i = 0;i < nx; i ++)
    y3 += dummy1[i];

  y3 -= (dummy1[0] + dummy1[nx - 1]) / 2;
  y3 *= dx;
  return y1 + y2 + y3; 
}

/*the inverse of fermi dirac integral for order 0 and 1/2, see 
 * nanomos */
double anti_dummy(double x, double order) {
  double y;
  double x2,x3,x4;

  x = fabs(x);
  if (x < 1e-14) x = 1e-14;

  if (fabs(order) < 1e-5) // then order = 0
    y = log(exp(x) - 1);
  else if (fabs(order - 0.5) < 1e-5) //order = 0.5
  {
    x2 = x * x;
    x3 = x2 * x;
    x4 = x3 * x;
    y=log(x) + 3.53553e-1 * x - 4.95009e-3 * x2 + 1.48386e-4 * x3 - 4.42563e-6 * x4 ; 
  } else {
    cout << "wrong order in anti_dummy : order = " << order << endl;
    exit(1);
  }
  return y;
}

/*when aprox = false, we use fermid_() to get better precision, else 
 * order 1/2 can be compute aproximately, see nanomos*/
double fermi_integral(double x, double order, bool aprox) {
  double relerr = 1e-8;
  int ierr;
  double fd;
  double exp_fac, nu, zeta, nu_prime, zeta_prime;

  if (!aprox){
    fermid_(&order, &x, &relerr, &fd, &ierr);  
    return fd;
  }

  if (fabs(order) < 1e-5) // fermi dirac integral of order 0
    fd = log(1 + exp(x));
  else if (fabs(order - 0.5) < 1e-5) // copied from nanomos
  {
    exp_fac = exp(-0.17 * pow(x + 1.0, 2));
    nu = pow(x,4) + 50.0 + 33.6 * x * (1.0 - 0.68 * exp_fac);
    zeta = 3.0 * sqrt(PI) / (4.0* pow(nu, 0.375));
    fd=exp(x) / (1.0 + zeta * exp(x));	
  } else if (fabs(order + 0.5) < 1e-5) {
    exp_fac = exp(-0.17 * pow(x + 1.0, 2));
    nu = pow(x,4) + 50.0 + 33.6 * x * (1.0 - 0.68 * exp_fac);
    zeta = 3.0 * sqrt(PI) / (4.0* pow(nu, 0.375));
    nu_prime = 4 * pow(x,3) + 33.6 - 22.848 * exp_fac * (1 - 0.34 * (x + x * x));
    zeta_prime = -(9 * pow(PI, 0.5) / 32) * pow(nu, -11.0/8) * nu_prime;
    fd = (exp(-x) - zeta_prime) / pow(exp(-x) + zeta, 2);
 } else {
    cout << "wrong order in fermi_integral(): order = " << order << endl;
    exit(1);
  }
  return fd;
}

void cp_to_buf(double * dbl_buf, int * int_buf, Particle par) {
  dbl_buf[0] = par.x;
  dbl_buf[1] = par.y;
  dbl_buf[2] = par.z;
  dbl_buf[3] = par.kx;
  dbl_buf[4] = par.ky;
  dbl_buf[5] = par.kz;
  dbl_buf[6] = par.charge;
  dbl_buf[7] = par.energy;
  dbl_buf[8] = par.left_time;
  
  int_buf[0] = par.i;
  int_buf[1] = par.j;
  int_buf[2] = par.k;
  int_buf[3] = par.itet;
  int_buf[4] = par.isym;
  int_buf[5] = par.par_type;
  int_buf[6] = par.flag;
  int_buf[7] = par.seed;
  int_buf[8] = par.par_id;
}

Particle make_particle(double * dbl_buf, int * int_buf) {
  Particle par;
  par.x  = dbl_buf[0];
  par.y  = dbl_buf[1];
  par.z  = dbl_buf[2];
  par.kx = dbl_buf[3];
  par.ky = dbl_buf[4];
  par.kz = dbl_buf[5];
  par.charge = dbl_buf[6];
  par.energy = dbl_buf[7];
  par.left_time = dbl_buf[8];
  par.i = int_buf[0];
  par.j = int_buf[1];
  par.k = int_buf[2];
  par.itet = int_buf[3];
  par.isym = int_buf[4];
  par.par_type = int_buf[5];
  par.flag = int_buf[6];
  par.seed = int_buf[7];
  par.par_id = int_buf[8];
  return par;
}

void init_int_val() {
  double max, min, bot, dx, dxx, x, s, y, xx;
  int i, j, n, N;
  max = 30;
  min = -30;
  bot = -50;
  dx = 1e-2;
  n = (int)((max - min) / dx);
  int_val.resize(n);
  N = 1000;

  for (i = 0;i < n;i ++){
    x = min + dx * i;
    dxx = (x - bot) / N;
    s = 0;
    for (j = 0;j <= N; j ++){
      xx = bot + j * dxx;
      y = log(1 + exp(xx));
      s += y;
    }
    s = s - log(1 + exp(bot)) * 0.5 - log(1 + exp(x)) * 0.5;
    s *= dxx;
    int_val[i] = s;
  }
}

double find_int_val(double x) {
  double max, min, bot, dx;
  int index;
  max = 30;
  min = -30;
  bot = -50;
  dx = 1e-2;
  if ((x < min) || (x > max)) 
    return 0;
  index = (int) ((x - min) / dx);
  return int_val[index];
}

/*use the 3th parameter to scale the output data, scale = 1 by default */
 void print_p_data(Epetra_Vector *p_data, string filename, double scale) {
   int i,j,k;
   MPI_Status status;

   if (rank != 0)
     MPI_Recv(&j, 1, MPI_INT, rank - 1, 99 , MPI_COMM_WORLD, &status);

   ofstream ofile;

   filename = "./data/" + filename;
   ofile.open(filename.c_str(), iostream::app);

   double *p_data_value;

   p_data->ExtractView(&p_data_value);

   for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
       for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++){
       ofile << i << ' ' << j << ' ' << k << ' ' << p_data_value[P_LINDEX_ONE(i,j,k)] * scale << endl;
     }

   ofile.close();

   if (rank != size - 1)
     MPI_Send(&j, 1, MPI_INT, rank + 1, 99, MPI_COMM_WORLD);

 }

 void print_c_data(Epetra_Vector *c_data, string filename, bool ghost, double scale) {
   int i,j,k;
   MPI_Status status;
   double *c_data_value;
   int jb, je;

 if (rank != 0)
     MPI_Recv(&j, 1, MPI_INT, rank - 1, 99 , MPI_COMM_WORLD, &status);

   ofstream ofile;

   filename = "./data/" + filename;

   ofile.open(filename.c_str(), iostream::app);

   c_data->ExtractView(&c_data_value);
   jb = c_jbegin;
   je = c_jend;
   if (ghost){
     jb = c_jbegin_ghost;
     je = c_jend_ghost;
   }

   for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
     for (j = jb; j <= je; j ++){
       ofile << i << ' ' << j << ' ' << k << ' ' << c_data_value[C_LINDEX_GHOST_ONE(i,j,k)] * scale << endl;
     }

   ofile.close();

   if (rank != size - 1)
     MPI_Send(&j, 1, MPI_INT, rank + 1, 99, MPI_COMM_WORLD);

 }


 void print_c_data(Epetra_IntVector *c_data, string filename, bool ghost, double scale) {
   int i,j,k;
   MPI_Status status;
   int *c_data_value;
   int jb, je;
   if (rank != 0)
     MPI_Recv(&j, 1, MPI_INT, rank - 1, 99 , MPI_COMM_WORLD, &status);

   ofstream ofile;

   filename = "./data/" + filename;

   ofile.open(filename.c_str(), iostream::app);

   c_data->ExtractView(&c_data_value);

   jb = c_jbegin;
   je = c_jend;
   if (ghost){
     jb = c_jbegin_ghost;
     je = c_jend_ghost;
   }

   for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
     for (j = jb; j <= je; j ++){
       ofile << i << ' ' << j << ' ' << k << ' ' << c_data_value[C_LINDEX_GHOST_ONE(i,j,k)] * scale << endl;
     }

   ofile.close();

   if (rank != size - 1)
     MPI_Send(&j, 1, MPI_INT, rank + 1, 99, MPI_COMM_WORLD);

 }

/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

double ran1(int * idum)
{
        int j,k;
	static int iv[NTAB],iy=0;
	void nrerror();
        static double NDIV = 1.0/(1.0+(IM-1.0)/NTAB);
        static double RNMX = (1.0-EPS);
        static double AM = (1.0/IM);

	if ((*idum <= 0) || (iy == 0)) {
		*idum = MAX(-*idum,*idum);
                for(j=NTAB+7;j>=0;j--) {
			k = *idum/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;
			if(*idum < 0) *idum += IM;
			if(j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = *idum/IQ;
	*idum = IA*(*idum-k*IQ)-IR*k;
	if(*idum<0) *idum += IM;
	j = (int) iy*NDIV;
	iy = iv[j];
	iv[j] = *idum;
	return MIN(AM*iy,RNMX);
}
#undef IA 
#undef IM 
#undef IQ
#undef IR
#undef NTAB
#undef EPS 
#undef MAX
#undef MIN

int ran_seed[1];

double Random() {
  return ran1(ran_seed);
}

void mysrand(int seed) {
  ran_seed[0] = -seed;
}

double random_fermi_fun(double Ef, double Emax) {
  double scale;
  double a = log(exp(-Ef) / (1 + exp(-Ef)));
  scale = log(1 - fermi_fun(Emax - Ef)) - a;
  return Ef - log(exp(-Random() * scale - a) - 1);
}

double fermi_fun(double x) {
  if (x > 40)
    return exp(-x);
  if (x < -40) 
    return 1;
  return 1.0 / (1 + exp(x));
}


