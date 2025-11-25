/**
 * @file mcmodel.h
 * @author Wei Zhang
 * @date 2010-4-13
 * @brief the head file of class MeshQuantities 
 * */

#ifndef __mcmodel__
#define __mcmodel__

#include "Constant_GlobalVariable_MathFunction.h"
#include "Band.h"
#include "myStatistic.h"
#include "utils.h"
#include "defs.h"

#include <libgen.h>

extern Band band;

class Amesos_PoissonSolver;
class PoissonSolver;

/* -------------------------------------------------------------------------- */
/** @brief core class for the simulation
 */
/* ---------------------------------------------------------------------------- */
class MeshQuantities {

  const Epetra_MpiComm * Comm; ///< used by trilinos communication

  PoissonSolver * poisson_solver;

  friend class PoissonSolver;

  Epetra_Vector * p_potential; 
  Epetra_Vector * c_field;
  Epetra_Vector * c_volume;
  
  Epetra_Vector * c_donor;
  Epetra_Vector * c_acceptor;
  Epetra_Vector * c_init_hole_charge;
  Epetra_Vector * c_init_electron_charge;
  Epetra_Vector * c_da;
  Epetra_Vector * p_dop;
  Epetra_Vector * p_volume;
  Epetra_Vector * p_vadd;
  Epetra_Vector * c_work1, *c_work2, *c_work3, * p_work1, *p_work2;
  
  Epetra_Vector * stat_vxE, * stat_vxH, *stat_vzE, *stat_vzH,*stat_vyE, *stat_vyH, *stat_energyE, *stat_energyH, *stat_chargeE, *stat_chargeH, *stat_pot;

  Epetra_Vector * stat_pot_saved, *stat_qc_pot_saved , * stat_chargeE_saved, * stat_chargeH_saved;

  Epetra_Vector * stat_qc_pot, *stat_ec, *stat_enum, *stat_e_heat, *stat_h_heat;
  
  Epetra_Vector * p_qc_range_x1, * p_qc_range_x2, * p_qc_range_y1, * p_qc_range_y2, * p_qc_range_z1, * p_qc_range_z2, * p_qc_Eb;

  vector<double> qclx, qcly, qclz;

  Epetra_Vector * p_charge_fac;

  Epetra_Import * p_Importer, * c_Importer, *p_Identity_Importer, *p_qc_Importer;

//  Epetra_Export * p_heat_Exporter;

  Epetra_IntVector * c_material;
  Epetra_IntVector * c_desired_electron_number;
  Epetra_IntVector * c_desired_hole_number;
  Epetra_IntVector * c_motion_rule;
  Epetra_IntVector * c_scatter_type;
  Epetra_IntVector * c_attached_contact;
  Epetra_IntVector * c_electron_num;
  Epetra_IntVector * c_hole_num;
  Epetra_IntVector * p_material;

  Epetra_Time * t_time;

  Epetra_Vector * p_pot, *p_tot_pot, * p_poisson_pot, *p_poisson_pot_saved, *p_rhs, * new_Ec;

  Epetra_Vector *c_field_x, *c_field_y, *c_field_z; 

  Epetra_Vector *c_h_field_x, *c_h_field_y, *c_h_field_z; 

  Epetra_Vector * p_nq;
  
  Epetra_Vector * p_fermi_level;

  Epetra_Vector * p_qc_fermi_level;

  Epetra_Vector * p_lnpot, * p_qcoeff, * p_quantum_stat_pot, * p_quantum_aver_qpot;

  Epetra_Vector * p_qc_pot;

  Epetra_IntVector * p_qc_mat;
  
  Epetra_Vector * c_par_charge;

  Epetra_Vector * p_par_charge, *old_p_par_charge;

  Epetra_Vector * c_electron_charge;

  Epetra_Vector * c_hole_charge;

  Epetra_IntVector * c_InSurfRegion, *c_quantumRegion;

  Epetra_IntVector * c_nearestSurf, *c_EeffDirection;

  Epetra_Vector * c_e_roughnessRate, * c_e_surfPhononRate;
  Epetra_Vector * c_h_roughnessRate, * c_h_surfPhononRate;

  Epetra_Vector * p_electron_heat, * p_hole_heat, *p_e_heat_weight, *p_h_heat_weight;

  int num_local_y;
  
  vector<list<Particle> >  par_list;

  list<Particle> * current_par_list;

  list<Particle>::iterator  par_iter;
  

  int local_par_num, par_num, inject_par_num;

  double sum_charge[2];
  double sum_charge_proc[2];
  
  
  string grid_file_name;
  string device_file_name;

  vector<user_cmd> cmd_list;
  vector<Contact> contact;

  vector<int> lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;

  vector<vector<vector<double> > >  subbands, density_sub;
  vector<vector<vector<vector<vector<double> > > > > subband_vec;

  myStatistic sttt;

  int step, total_step, electron_number, hole_number, default_electron_num, default_hole_num;

  /**
   * @brief device temperature
   * 
   */
  double device_temperature;

  bool Flag_MultipleRefresh;

  bool Flag_compute_potential;

  bool Flag_compute_heat, flag_heat;

  bool Flag_QuantumCorrection;

  bool Flag_LaterQC;

  bool Flag_SPE;

  bool Flag_SurfaceScatter, Flag_SurfaceRoughnessScatter, Flag_SurfacePhononScatter, Flag_SurfaceImpurityScatter;

  bool Flag_SelfScatter;

  bool Flag_Catch;

  bool Flag_GetTetTime;

  bool Flag_NonLinearPoisson;
  
  bool Flag_restart;
  
  bool Flag_calSurfscatt;
  
  bool Flag_test;

  int restart_step;

  int no_QC_step;

  int stat_step, stat_heat_step, qc_stat_step, debug_print_step, quantum_start_step,quantum_stat_step, quantum_print_steps, quantum_step , quantum_aver_step;

  int max_subband, valley_num;

  int mr_step, heat_steps;

  string output_dir, bs_path, restart_filename;
  string kloem_table_file, kloab_table_file, ktoem_table_file, ktoab_table_file,
         klaem_table_file, klaab_table_file;

  double Rpar, Rsqr, pckill;

  double dt;

  double inject_const;

  double frickel, Tn, QuantumPotentialCoef;
  
  double SurfSc_ail, SurfSc_delta, SurfSc_XIph, SurfSc_theta, SurfSc_gama, SurfSc_Nbmod, SurfSc_Npf, SurfSc_Pft,
    SurfSc_Pfn, SurfSc_Kpha, SurfSc_Kt, SurfSc_Fsimp, SurfSc_Rshmin, SurfSc_Rshmax, SurfSc_GAMMAn;
  
  int SurfaceType[MNSurface];
  double SurfacePosition[MNSurface];
  int SurfaceDir[MNSurface];

	
  double energy, kx,ky,kz, charge;
  int iband, itet, isym, itetdir;
  double vx,vy,vz, left_time;
  int par_type;
  int icell, jcell, kcell;
  double x, y,z;
  int idir;
  int seed;
  int par_id;

  double NumSurface;

  double Ex, Ey, Ez, Rho, DA;

  /**
   * @param TetTf:    Tetra Time
   * @param CellTf:   Time of Flight in Single Cell
   * @param PhScTf:   Phonon Scatter time
   * @param ImpScTf:  Impact Scatter Time
   * @param SurfScTf: Surface Scatter Time
   */
  double TetTf, CellTf, PhScTf, ImpScTf, SurfScTf, Tf;

  double ImpScGamma, SurfScGamma;

  int catch_par_num, unfinish_par, gen_par, mr_gen_num;

  double SurfRoughnessScRate, SurfPhononScRate;

  double * lcurrent, * rcurrent;

  void get_state();
  
  void getK();
  
  void out_wedge(double &xout,double &yout,double &zout,
                double xin,double yin,double zin);

  void in_wedge(double &xout,double &yout,double &zout,
                double xin,double yin,double zin);
  
  void GetV();

  void GetEnergy();

  void GetisymBZ();

  void Getisym();

  void TALIND(int);

  void ESCATII(int);

  void EPSCAT(int);

  void EFPSCAT(int);

  double FEKLOAB(double);

  double FEKLOEM(double);

  double FEKLAAB(double);

  double FEKLAEM(double);
  
  double FEKTOAB(double);

  double FEKTOEM(double);

  double EOLINT(double);

  void  FPFSAB();

  void FPFSEM();

  void init_deep();

  void ElectronPhononScatter();

  void ElectronImpurityScatter();

  void EleImpScUpdateK();

  void HolePhononScatter();

  void HSCATII(int);

  void HPSCAT(int);
  
  void ParticleSurfaceScatter();

  void HoleSurfacePhononScatter();

  void ElectronSurfacePhononScatter();

  void GetSurfRoughnessPhononScRate();

  void HitTet();

  int HitCell();

  void Generate();

  void Period();

  void Diffuse();

  void Reflect();

  void CatchAtGate();

  void CatchAtContact();

  void Pass();
  
  double TetTime();

  double CellTime();

  double GetImpScGamma();

  double GetImpScRate(double);

  void GetSurfScRate();
  
  void OutPar(list<Particle>::iterator );

  void OutPar(Particle *);

  void InPar(list<Particle>::iterator );
    
  double kloem[210][130];
  double kloab[210][120];
  double ktoem[60][170];
  double ktoab[60][160];
  double klaem[60][170];
  double klaab[60][160];
  
  void read_grid_file();

  void read_device_file();

  void get_cube_range(ifstream &ifile, double * pos, int * index);

  int distance2index(double pos, vector<double> &cut);

  void scaling();

  void init_phpysical_parameter(char * filename);

  void getInputData(char *);

  void init_wedge();

  void init_epetra_map_vector();
  
  void init_vector_entry(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK, int EndK,
                         double entry_value, Epetra_Vector * c_vector);

  void init_vector_entry(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK, int EndK,
                         int entry_value, Epetra_IntVector * c_vector);
  
  void init_by_user_input();

  void init_cell_data();

  void init_surface_roughness();

  void init_point_data();

  void init_particle_data();

  void compute_total_par_num();

  void compute_par_num();

  void print_particle_data(string);

  void print_subband_energy() ;

  void init_poisson_matrix();

  void init_nonlinear_poisson();
  
  void SetMRCube(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK, int EndK,
                 int mr0,int mr1,int mr2,int mr3, int mr4, int mr5);
  
  void SetMRPlane(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK, int EndK,

                  int idn,int mr);
  
  void compute_field();

  void compute_fermi_level(bool, Epetra_Vector *);

  void compute_cell_charge();

  void nonlinear_poisson_solver(Epetra_Vector *, bool);

  void linear_poisson_solver() ;

  void density();
  
  double Ec_diff();
  
  void find_maxloc(double &, int &, int , int);

  void quantumPotential();

  void save_poisson_pot();

  void init_effective_potential();

  double effective_potential(int pi, int pj, int pk);

  void init_p_mat();

  void init_ep_range();

  void check_par_number();
  
  void pot2field(int, Epetra_Vector *);

  void pot2field_for_hole(Epetra_Vector *);
  
  void update_particle();

  void dump_par_info();

  void dump_par_info(Particle &);

  void trace_back();

  void dump_time_info();

  void dump_cell_info(int i, int j,int k);

  void particle_fly();

  void migrate_particle();

  void local_migrate();

  void global_migrate();

  void reduce_add(Epetra_Vector *, int );
  
  void clear_ghost_par();
  
  void compute_rhs();

  void statistic();

  void set_stat_zero();

  void fill_ghost_cell();

  void output_stat();

  void current_scatter_info();

  void inject_particle(bool source, int *, double *);

  void select_kstate(Particle * iter, int dir);

  void select_kstate_fermi_dirac(Particle * iter, int dir, double );

  int inject_cell_num(int i, int j, int k,  double * average_charge);

  void adjust_sd_charge();

  void adjust_charge(bool);

  void particle_to_density(Epetra_Vector *, int);

  int carrier_inject();

  Epetra_CrsMatrix * init_matrix();

  void print_qc_pot() ;

  void quantum_stat_pot() ;

  void clear_quantum_stat_pot();

  void MultipleRefresh();

  bool gen_new_par_mr(vector<Particle> & , int , vector<Particle> & );
  
  void read_potential();

  void read_par_info_for_restart(Epetra_Vector *, Epetra_IntVector * , Epetra_Vector * , Epetra_IntVector * );

  void output_for_restart();

  void stat_heat_single_step();

  void compute_heat();

  void heat_field();

  void heat_density();

  void heat_to_point(double, Epetra_Vector *);

  void read_device_input_temperature(char *);

public:
  
  MeshQuantities();

  void initialize(char *);

  void run();

};

#endif

