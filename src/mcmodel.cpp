/**
 * @file mcmodel.cpp
 * @author Wei Zhang
 * @date 2010-4-13
 * @brief the major file of class MeshQuantities 
 * */

#include "mcmodel.h"
#include "PoissonSolver.h"

#include <sys/stat.h>
#include <cerrno>
#include <cstring>
#include <unistd.h>

/* -------------------------------------------------------------------------- */
/** @brief 能带结构的对象
 */
/* ---------------------------------------------------------------------------- */

Band band;

/** 
 * @brief 默认构造函数s
 */
MeshQuantities::MeshQuantities() {
  
  Comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  mpi_rank = Comm->MyPID();

  mpi_size = Comm->NumProc();

  t_time = new Epetra_Time(*Comm);
}

/* -------------------------------------------------------------------------- */
/** @brief 模拟之前初始化
 * 
 * @param FileName 用户输入的文件名
 * 
 * @return   
 */
/* ---------------------------------------------------------------------------- */
void MeshQuantities::initialize(char * FileName){

  /* 打印开始初始化，读取参数*/
  cout << "Initializing ... Reading input parameters from file: " << FileName << endl;
  
  /*初始化一些物理参数*/
  // 初始化参数，并对部分物理参数做了去单位化的处理
  init_phpysical_parameter(FileName);
  

  /* 读入用户提供的，仿真所使用的模拟参数 */
  getInputData(FileName);
  

  /* 读入物理模型的网格数据 */
  // 读取 grid 文件 (lgrid.txt) 中的格点信息
  read_grid_file();
  

  /*define the vectors*/
  init_epetra_map_vector();    //暂时不管

  

  /*load used-provided file*/
  /** 
   * @brief 读取输入的器件的具体信息 (ldg.txt) 
   * ldg.txt 中的每一行命令会被存储为一个 user_cmd 对象 cmd 
   * 然后会将所有 user_cmd 对象存储在一个 cmd_list 之中
   */
  read_device_file();

  // 某个能量区间内，不同的 k 空间采样点区间对应的概率密度
  init_deep();

  /**
   * @brief init the vectors according to device file
   * 从 cmd_list 之中，逐行取出存储的 user_cmd 对象，即输入文件中的各个命令
   * 然后根据命令的关键字对应的编号值，为程序中的 vector 进行相应的赋值
   */
  init_by_user_input();

  /* 去单位化 */
  scaling();

  

  /*init band structure*/
  /**
   * @brief init band structure
   * 为很多参数做了去单位化
   * 读取了 band structure 文件中的内容
   * 计算了态密度
   * 计算了电子、空穴和氧化层电子 的 声子散射+电离杂质散射 的散射率中除态密度之外的部分
   * 计算了声子散射+电离杂质散射 的散射率之和，是直接求和。
   * 计算了本征载流子密度
   */
  band.IELEC(bs_path);

  /*init cell data*/
  init_cell_data();

  if (Flag_SurfaceScatter)
    init_surface_roughness();

  /*init point data*/
  init_point_data();
  
  init_particle_data(); /*init particles according to nuetrality condition*/
  
  /*init poisson matrix*/
  init_poisson_matrix();

  

  // /*initialize for quantum correction*/
  // if (Flag_QuantumCorrection){
  //   clear_quantum_stat_pot();
  //   quantum_step = 0;
  //   quantum_aver_step = 0;
  //   stat_qc_pot->PutScalar(0);
  // }

  rcurrent = (double *) malloc(p_numy * sizeof(double));
  lcurrent = (double *) malloc(p_numy * sizeof(double));

  for (int i = 0;i < p_numy; i ++) {
    rcurrent[i] = 0;
    lcurrent[i] = 0;
  }

  init_effective_potential();

  init_p_mat();

  init_ep_range();

   if (!Flag_compute_potential)
     read_potential();

  
}
/* --------------------------------------------------------------------- */
/**
 * @brief read device temperature
 */
void MeshQuantities::read_device_input_temperature(char *filename) {

    Trilinos_Util::InputFileReader fileReader(filename);
    fileReader.ReadFile();

    device_temperature = fileReader.Get("Temperature", 300.0);
}

/* -------------------------------------------------------------------------- */
/** @brief 调试用， 不求解 poisson 方程，而从文件读入电势.
 * 
 */
/* ---------------------------------------------------------------------------- */
void MeshQuantities::read_potential() {
  int i,j,k;
  double tmp;
  double * pot;
  FILE * fin;
  int flag;

  fin = fopen("./input_pot", "r");

  p_poisson_pot->ExtractView(&pot);

  do {
    flag = fscanf(fin, "%d %d %d", &i,&j,&k);
    if (flag == 3) {
      fscanf(fin, "%lf", &tmp);
      if ((j >= p_jbegin_nonoverlap) && (j <= p_jend_nonoverlap))
         pot[P_LINDEX_ONE(i,j,k)] = tmp / pot0;
    }
    else break;
  } while(1);
  fclose(fin);
}
/* -------------------------------------------------------------------------- */
/** @brief 调试用，从源漏端注入粒子， 废除.
 * 
 */
/* ---------------------------------------------------------------------------- */

int MeshQuantities::carrier_inject() {

  int tot_inject = 0;
  
  int inject_num[2], tot_inject_num, drain_inject_num[2], source_inject_num[2];

  double inject_charge[2], source_inject_charge[2];

  MPI_Status status;

    /* inject particles from source contact */
  if (p_jbegin == 0){ 
    /*true means source injection, number and charge are returned*/
     inject_particle(true, source_inject_num, source_inject_charge);
/*
     contact[0].NumParGen += source_inject_num[0] + source_inject_num[1];
     contact[0].CharGen += source_inject_charge[0] + source_inject_charge[1];
     */
  }       

    /* inject particles from drain contact */
  if (p_jend == p_numy - 1){

     /* false means drain injection, number and charge are returned */
     inject_particle(false, inject_num, inject_charge);
/*
     contact[1].NumParGen += inject_num[0] + inject_num[1];
     contact[1].CharGen += inject_charge[0] + inject_charge[1];
     */

     /*for debug */
     // cout << inject_charge[0] << ' ' << inject_charge[1] << endl;

     MPI_Send(inject_num, 2, MPI_INT, 0, 99, MPI_COMM_WORLD);
  }

   if (mpi_rank == 0) {

     cout << "carrier inject:\n ";

     tot_inject_num = 0;

     MPI_Recv(drain_inject_num, 2, MPI_INT, mpi_size - 1, 99 , MPI_COMM_WORLD, &status);

     /* index 0 denotes electron number , while 1 denotes
      * hole number */
     tot_inject_num += source_inject_num[0] + source_inject_num[1];
 
     /* index 0 denotes electron number , while 1 denotes
      * hole number */
     tot_inject_num += drain_inject_num[0] + drain_inject_num[1];

     cout << " source : " << source_inject_num[PELEC] << " electrons, " 
       << source_inject_num[PHOLE] << " holes\n";

      /*report */
     cout << "  drain : " << drain_inject_num[PELEC] << " electrons, "
	 << drain_inject_num[PHOLE] << " holes\n";
   }

   /* have a rest here*/
   MPI_Barrier(MPI_COMM_WORLD);

   return tot_inject_num;
}

/* -------------------------------------------------------------------------- */
/** @brief 调试用，找电势最大值的位置, 废除.
 * depracated
 * 
 * @param Ec_peak 
 * @param j_peak 
 * @param val 
 * @param band 
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::find_maxloc(double  & Ec_peak, int & j_peak, int val, int band) {

  struct {
    double value;
    int index;
  } in, out;

  int j;

  in.value = -100;
  for (j = p_jbegin_nonoverlap;j <= p_jend_nonoverlap; j ++)
    if (in.value < subbands[P_SHIFT_NONOVERLAP_Y(j)][val][band]) {
      in.value = subbands[P_SHIFT_NONOVERLAP_Y(j)][val][band];
      in.index = j;
    }

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); 

  Ec_peak = out.value;
  j_peak = out.index;
}
/* -------------------------------------------------------------------------- */
/** @brief 计算 fermi 能级, 废除.
 * 
 * @param pot_based 
 * @param p_pot_vec 
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::compute_fermi_level(bool pot_based, Epetra_Vector * p_pot_vec){

  vector<double> Fn;
  int i,j,k;
  double Fn_slope = (Vs - Vd) / (ly[p_gend] - ly[p_gbegin]);
  double *fermi_val;
  double *pot;

  p_fermi_level->ExtractView(&fermi_val);

  if (!pot_based) {
    Fn.resize(p_num_local_nonoverlap_y);

    for (j = p_jbegin_nonoverlap;j <= p_jend_nonoverlap; j ++)
      if (j < p_gbegin) 
	Fn[j - p_jbegin_nonoverlap] = - Vs ;
      else if (j > p_gend)
	Fn[j - p_jbegin_nonoverlap] = - Vd ; 
      else 
	Fn[j - p_jbegin_nonoverlap] = - Vs + (ly[j] - ly[p_gbegin]) * Fn_slope;

    for (i = p_tox; i <= p_box; i ++)
      for (j = p_jbegin_nonoverlap;j <= p_jend_nonoverlap; j ++)
	fermi_val[P_LINDEX_ONE(i,j,k)] = Fn[j - p_jbegin_nonoverlap] / pot0;
  } else {
    p_pot_vec->ExtractView(&pot);
    double *p_par_charge_value;
    p_par_charge->ExtractView(&p_par_charge_value);
    for (i = p_tox; i <= p_box; i ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      if (NOT_GHOST_BC(j)){
	  fermi_val[P_LINDEX_ONE(i,j,k)]
	    = anti_dummy(p_par_charge_value[P_LINDEX_ONE(i,j,k)] * conc0 / Nc, fermi_order) + pot[P_LINDEX_ONE(i,j,k)];
      }
  }
/*
#ifdef DEBUG_DATA_PRINT
    string filename;
    filename = getFileName("fermi_level", step);
    print_p_data(p_fermi_level,filename, pot0);
#endif
*/
}
/* -------------------------------------------------------------------------- */
/** @brief 根据粒子信息统计浓度分布.
 * 
 * @param p_par 保存计算出的浓度.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::particle_to_density(Epetra_Vector * p_par, int par_flag) {

    int i,j,k;

    double ratio_x, ratio_y,ratio_z;

    double *vol, *pot, cc;

    double *p_par_charge_value;

    list<Particle> * c_par_list;

    list<Particle>::iterator iter;

    p_volume->ExtractView(&vol);

    p_par->PutScalar(0);

    p_par->ExtractView(&p_par_charge_value);
 
     /* loop for each cell */ 
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	for (j = c_jbegin_ghost; j <= c_jend; j ++)
//        if (NOT_GHOST_CELL(j))
	{
	  /* get the cell's particle list */	
	  c_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];

	  /*for each particle */
	  for (iter = c_par_list->begin(); iter != c_par_list->end(); iter ++) {

            if ((par_flag != 2) && (par_flag != iter->par_type))
              continue;

	    cc = iter->charge;

	    ratio_x = (iter->x - lx[i]) / dx[i];
	    ratio_y = (iter->y - ly[j]) / dy[j];
	    ratio_z = (iter->z - lz[k]) / dz[k];

	    if (!(BETWEEN01(ratio_x) && BETWEEN01(ratio_y) && BETWEEN01(ratio_z))) {
	      err_message(WRONG_CELL, "particle to density");
              dump_par_info(*iter);
	    }
	    /* distribute the contribution to cell's each node */
	    if (j >= p_jbegin) {
	      p_par_charge_value[P_LINDEX_ONE(i,j,k)] += cc * (1 - ratio_y) * (1 - ratio_x) * (1 - ratio_z);
	      p_par_charge_value[P_LINDEX_ONE(i + 1, j,k)] += cc * (1 - ratio_y) * ratio_x * (1 - ratio_z);
	      p_par_charge_value[P_LINDEX_ONE(i,j,k+1)] += cc * (1 - ratio_y) * (1 - ratio_x) * ratio_z;
	      p_par_charge_value[P_LINDEX_ONE(i + 1, j,k + 1)] += cc * (1 - ratio_y) * ratio_x * ratio_z;
	    }
	    if (j + 1 <= p_jend_nonoverlap) {
	      p_par_charge_value[P_LINDEX_ONE(i,j + 1,k)] += cc * ratio_y * (1 - ratio_x)* (1 - ratio_z);
	      p_par_charge_value[P_LINDEX_ONE(i + 1, j + 1,k)] += cc * ratio_y * ratio_x* (1 - ratio_z);
	      p_par_charge_value[P_LINDEX_ONE(i,j + 1,k+1)] += cc * ratio_y * (1 - ratio_x) * ratio_z;
	      p_par_charge_value[P_LINDEX_ONE(i + 1, j + 1,k+1)] += cc * ratio_y * ratio_x * ratio_z;

	    }
	  }
	}
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap;j <= p_jend_nonoverlap; j ++)
//	if (NOT_GHOST_BC(j))
      {
      if (fabs(vol[P_LINDEX_ONE(i,j,k)]) > 0)
	  p_par_charge_value[P_LINDEX_ONE(i,j,k)] = p_par_charge_value[P_LINDEX_ONE(i,j,k)] / vol[P_LINDEX_ONE(i,j,k)]; 
      }
}

/* -------------------------------------------------------------------------- */
/** @brief 主要的模拟函数，包括求解 poisson 方程和模拟粒子两部分
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::density() {
  double *fermi_val, *charge_fac, *qc_fermi_val;
  double *p_par_charge_value;
  double *pot;

  // charge_fac is pointing to p_charge_fac
  // 通过指针获取电荷系数数组
  p_charge_fac->ExtractView(&charge_fac);
  // 将 p_par_charge 数组清零
  p_par_charge->PutScalar(0);
  // 获取粒子所带电荷值的数组的指针
  p_par_charge->ExtractView(&p_par_charge_value);
 
  /**
  * @brief 填充 ghost cell 区域的数据，确保不同的 MPI 进程间的边界数据能够正确交换
  */
  fill_ghost_cell();

  // 计算和统计上一步的结果
  /**
   * @brief 统计上一步的模拟结果，计算各种物理量的统计值
   * 物理量包括：速度分量、电荷量、载流子能量
   */
  statistic();

  // 计算每个粒子的单位电荷
  //  粒子的体积为其所属的各个单元的体积之和 * 1/8
  /**
   * @brief 
   * 根据每个 cell 中的粒子的电荷量分配到每个 cell 的各个顶点上，
   * 分配时按照空间权重分配
   */
  particle_to_density(p_par_charge, 2);

  /**
  * @brief 基于求解得到的电势与电荷分布计算电场
  */
  compute_field();

  compute_cell_charge();

//  inject_par_num = carrier_inject();

//  particle_to_density(p_par_charge);

  clear_ghost_par();

    //adjust_sd_charge();

  if (Flag_SurfaceScatter)
    GetSurfRoughnessPhononScRate();

  t_time->ResetStartTime();

  update_particle();

  if (mpi_rank == 0)
    cout << "time used: " << t_time->ElapsedTime() << "s" << endl;

#ifdef DEBUG_DATA_PRINT
  /*
  if (step % debug_print_step == 0) {
    string filename;
    filename = getFileName("density", step);
    print_p_data(p_par_charge,filename, conc0);
  }
  */
#endif
}
/* -------------------------------------------------------------------------- */
/** @brief 计算相邻两步电势的差, 废除. 
 * 
 */
/* ---------------------------------------------------------------------------- */

double  MeshQuantities::Ec_diff() {
  Epetra_Vector tmp(*new_Ec);
  double tmp_norm;
  int i,j,k;
  double * tmp_value, * old_Ec_val, * new_Ec_val;

  tmp.ExtractView(&tmp_value);
  p_poisson_pot->ExtractView(&old_Ec_val);
  new_Ec->ExtractView(&new_Ec_val);

  for (i = p_ibegin; i <= p_iend; i ++)
    for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      tmp_value[P_LINDEX_ONE(i,j,k)] = old_Ec_val[P_LINDEX_ONE(i,j,k)] - new_Ec_val[P_LINDEX_ONE(i,j,k)];
  tmp.Norm2(&tmp_norm);
  return tmp_norm;
}
/* -------------------------------------------------------------------------- */
/** @brief 主要的模拟过程. 
 *  @return
 */
/* ---------------------------------------------------------------------------- */
void MeshQuantities::run() {

 //cout << "Running" << endl;
  if (mpi_rank == 0) {
    cout << "Entering run(): total_step = " << total_step
         << ", restart step = " << restart_step << endl;
  }

  if (Flag_restart) // run from the previous simulation result
    step = restart_step;
  else 
    step = 0;

  /*output for integration*/
  print_p_data(p_volume,"pvolume", conc0);

  flag_heat = false;

  /* simulation loop */
  for (; step < total_step; step ++) {

    if (mpi_rank == 0) {
      cout << "************ step = " << step << " ********" << endl;
      cout << par_num << " particles"  << endl;
    }

    /* self-consistent iteration for one step */
    density();

    /* this should be impossible , but still check it */
    if (par_num == 0) {
      cout << "par_num become zero" << endl;
      exit(0);
    }

    /* deal with Multiple Refresh */
    mr_gen_num = 0;

    if ((Flag_MultipleRefresh) && (step % mr_step == mr_step - 1)){

      MultipleRefresh();

      int tot_mr_gen_num;

      MPI_Allreduce(&mr_gen_num, &tot_mr_gen_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if (mpi_rank == 0)
        cout << "MultipleRefresh Gen : " << tot_mr_gen_num << endl;
    }

    /*particle should be conserved, check it */
    // Debug helper: audit particle distribution across valid/ghost cells
    if (mpi_rank == 0) Audit_Particles();
   check_par_number();

   /* output statistical result, every stat_step, this may be wrong
    * for the first statistical result, if run from restart */
    if (step % stat_step == (stat_step - 1)){
      output_stat();
      /* output restart result*/
      output_for_restart();
      /*set to zero for the following statistics*/
      set_stat_zero();
    }
  }

  if (Flag_compute_heat) 
    compute_heat();

}

/* -------------------------------------------------------------------------- */
/** @brief 从源漏端注入粒子, 废除.
 * 
 * @return the number of particles to be injected for cell (i,j) with type iptype, with average
 * charge equals to average_charge * 
 */
/* ---------------------------------------------------------------------------- */

int MeshQuantities::inject_cell_num(int i, int j,int k, double * aver) {

   double * par_charge, *ccharge_val,*c_da_value;
   double *volume_value;
   double tot_charge;
   int old_par_num ;
   int  *par_num;
   int num = 0;
   int maxGen = 100;
   double da_vol;

   c_electron_charge->ExtractView(&ccharge_val);
   old_par_num = par_list[C_LINDEX_GHOST_ONE(i,j,k)].size();
   c_da->ExtractView(&c_da_value);
   c_volume->ExtractView(&volume_value);

   da_vol = - c_da_value[C_LINDEX_GHOST_ONE(i,j,k)] * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];

   if (ccharge_val[C_LINDEX_GHOST_ONE(i,j,k)] > da_vol){

     tot_charge = da_vol - ccharge_val[C_LINDEX_GHOST_ONE(i,j,k)]; 

//   tot_charge = - band.Ef_to_cross_number(Ef - Ec) * dx[i] * dt; 

    *aver = ccharge_val[C_LINDEX_GHOST_ONE(i,j,k)] / old_par_num;

//     num = (int)round(tot_charge / *aver);

     num = (int) (tot_charge / *aver);

     if (num == 0) 
       return 0;
   
     if (num > maxGen) num = maxGen;

     *aver = tot_charge  / num;

   } else num = 0;

     return num;

  }


/* -------------------------------------------------------------------------- */
/** @brief 注入粒子.  
 * 
 * @param source true 表示源端注入, false 表示漏端注入.
 * @param total_inject_num 要注入的粒子数.
 * @param total_charge 要注入的总电荷数.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::inject_particle(bool source, int * total_inject_num, double * total_charge) {
  int i = 0,k;
  int j, cj;/*cell index*/
  double * pot, Ec, par_charge;
  int shift;
  double dope, Ej;
  int N_inject, ipar;
  double lamda;
  Particle newpar;
  double inject_charge;
  int dir;
  int * material;
  
  total_inject_num[0] = total_inject_num[1] = 0;
  total_charge[0] = total_charge[1] = 0;

  p_pot->ExtractView(&pot);
  c_material->ExtractView(&material);

  if (source) {
    /*source injection*/
    cj = 0;
    j = 0;
    Ej = fermi[0];
    dir = 1;
  }
  else {
    /*drain injection*/
    cj = c_numy - 1;
    j = p_numy - 2;
    Ej = fermi[1];
    dir = -1;
  }
  
  /*injection loop for each row, for silicon region */

    /*charge should be injected*/
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
        if (material[C_LINDEX_GHOST_ONE(i,cj,k)] == SILICON) {
    /* number of particles to be injected */
//	N_inject = inject_cell_num(i, cj,k, Ej, Ec, &par_charge);

	N_inject = inject_cell_num(i,cj,k, &par_charge);

	/*update based on the return results*/
	total_inject_num[0] += N_inject; 
	inject_charge = N_inject * par_charge;
	total_charge[0] += inject_charge;

	/*initialize for each injected particles */
	for (ipar = 0; ipar < N_inject; ipar ++){

	   seed = i * 100000 + cj * 1000 + ipar;
	   newpar.par_id = seed;
	   lamda = Random();
	   //             lamda = 0.5;
	   newpar.x = lamda * lx[i] + (1 - lamda) * lx[i + 1];

	   lamda = Random();
	   newpar.y = lamda * ly[j] + (1 - lamda) * ly[j + 1];

	   lamda = Random();
	   newpar.z = lamda * lz[k] + (1 - lamda) * lz[k + 1];

	   newpar.charge = par_charge;
	   newpar.i = i;
	   newpar.j = cj;
	   newpar.k = k;
	   newpar.par_type = 0;
	   newpar.seed = seed;
	   newpar.left_time = dt;

	   /*select kstate */
	   select_kstate(&newpar, 0);
//	   select_kstate_fermi_dirac(&newpar, dir, Ej - Ec);

	   /*add to list of cell particles*/
	   par_list[C_LINDEX_GHOST_ONE(i,cj,k)].push_back(newpar);
	}
    }
  }
/* -------------------------------------------------------------------------- */
/** @brief 统计数组清零
 * */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::set_stat_zero() {
  int i;
  stat_pot->PutScalar(0);
  stat_qc_pot->PutScalar(0);
  stat_vxE->PutScalar(0);
  stat_vxH->PutScalar(0);
  stat_vyE->PutScalar(0);
  stat_vyH->PutScalar(0);
  stat_vzE->PutScalar(0);
  stat_vzH->PutScalar(0);
  stat_chargeE->PutScalar(0);
  stat_chargeH->PutScalar(0);
  stat_energyE->PutScalar(0);
  stat_energyH->PutScalar(0);
  stat_ec->PutScalar(0);
  stat_enum->PutScalar(0);
  stat_e_heat->PutScalar(0);
  stat_h_heat->PutScalar(0);

  sttt.reset();
  for (i = 0;i < contact.size();i ++)
    contact[i].reset();

  for (i = 0;i < p_numy ; i ++)
    lcurrent[i] = 0;

  for (i = 0;i < p_numy; i ++)
    rcurrent[i] = 0;

  if (flag_heat)
    {
      p_electron_heat->PutScalar(0);
      p_hole_heat->PutScalar(0);
    }
}
/* -------------------------------------------------------------------------- */
/** @brief 计算并输出散射信息和电流值. 
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::current_scatter_info() {

  double gen_tmp[MNContact], catch_tmp[MNContact], energy_gen_tmp[MNContact], energy_catch_tmp[MNContact];
  int num_gen_tmp[MNContact], num_catch_tmp[MNContact];
  double total_gen[MNContact], total_catch[MNContact];
  double total_energy_gen[MNContact], total_energy_catch[MNContact];
  int total_genNum[MNContact], total_catchNum[MNContact];
  double * lcurrent_reduced, * rcurrent_reduced;

  int icont;
  int NumContact = contact.size();

  for (icont = 0 ; icont < NumContact; ++ icont) {
    gen_tmp[icont] = contact[icont].CharGen;
    energy_gen_tmp[icont] = contact[icont].EnergyGen;
    catch_tmp[icont] = contact[icont].CharCatch;
    energy_catch_tmp[icont] = contact[icont].EnergyCatch;
    num_gen_tmp[icont] = contact[icont].NumParGen;
    num_catch_tmp[icont] = contact[icont].NumParCatch;
  }

  MPI_Allreduce(gen_tmp, total_gen, NumContact, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(catch_tmp, total_catch, NumContact, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(energy_gen_tmp, total_energy_gen, NumContact, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(energy_catch_tmp, total_energy_catch, NumContact, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(num_gen_tmp, total_genNum, NumContact , MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(num_catch_tmp, total_catchNum, NumContact, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int phonon_scatter_tmp = sttt.phononScatter;
  int impurity_scatter_tmp = sttt.impurityScatter;
  int phonon_scatter, impurity_scatter;

  MPI_Allreduce(&phonon_scatter_tmp, &phonon_scatter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&impurity_scatter_tmp, &impurity_scatter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    ofstream ofile;
    string filename;
    filename = "./data/current";
    if (step == 0) 
      ofile.open(filename.c_str(),iostream::trunc); 
    else 
      ofile.open(filename.c_str(),iostream::app);

    ofile << "step = " << step << endl;
    
    for (icont = 0 ; icont < NumContact; ++ icont) {
      ofile << "icont = " << icont << endl;
      ofile << "CharGen = " << total_gen[icont] << ' '
          << "NumParGen = " << total_genNum[icont] << ' '
          << "CharCatch = " << total_catch[icont]  << ' '
          << "NumParCatch = " << total_catchNum[icont]  << endl
          << "EnergyGen= " << total_energy_gen[icont]  << endl
          << "EnergyCatch= " << total_energy_catch[icont]  << endl
          << "EnergyCurrent= " << (total_energy_catch[icont] - total_energy_gen[icont]) / (dt *  stat_step) * pot0  << endl
            << "Current = " << (total_gen[icont] - total_catch[icont]) /  (dt * stat_step)  * curr0 * 1e-2 << endl;
    }
    ofile << "times of phonon scatter : " << phonon_scatter << endl
          << "times of impurity scatter : " << impurity_scatter << endl;
    ofile.close();
  }

  rcurrent_reduced = (double *) malloc(p_numy * sizeof(double));
  lcurrent_reduced = (double *) malloc(p_numy * sizeof(double));

  MPI_Allreduce(rcurrent, rcurrent_reduced, p_numy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(lcurrent, lcurrent_reduced, p_numy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    ofstream ofile;
    string filename;
    filename = getFileName("./data/section_current", step);
    ofile.open(filename.c_str(),iostream::trunc); 
    int i;

    for (i = 0; i < p_numy; i ++)
      ofile << i << ' ' 
	    << rcurrent_reduced[i] / (dt * stat_step) * curr0 << ' ' 
	    << lcurrent_reduced[i] / (dt * stat_step) * curr0 << ' '  
	    << (rcurrent_reduced[i] - lcurrent_reduced[i]) / (dt * stat_step) * curr0 << endl;  
    ofile.close();
  }

  free(rcurrent_reduced);
  free(lcurrent_reduced);

}
/* -------------------------------------------------------------------------- */
/** @brief 从之前的结果开始计算.
 * 
 * @param e_charge 输入的浓度分布
 * @param e_num 每个单元格中电子数目.
 * @param h_charge 输入的每个单元格中空穴的浓度.
 * @param h_num 每个单元格中空穴的数目.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::read_par_info_for_restart(Epetra_Vector * e_charge, Epetra_IntVector * e_num, Epetra_Vector * h_charge, Epetra_IntVector * h_num) {

  ifstream restart_file;
  int n1, n2;
  int i,j,k;
  double c1, c2;
  int * en, * hn;
  double * ec, * hc;

  restart_file.open(restart_filename.c_str());

  e_charge->ExtractView(&ec);
  e_num->ExtractView(&en);

  h_charge->ExtractView(&hc);
  h_num->ExtractView(&hn);

  while (1) {

    restart_file >> i >> j >> k;

    if (i < 0) break;

    restart_file >> c1 >> n1 >> c2 >> n2;

    if ((i >= c_ibegin) && (i <= c_iend) && (j >= c_jbegin) && (j <= c_jend) && (k >= c_kbegin) && (k <= c_kend))
    {
      ec[C_LINDEX_GHOST_ONE(i,j,k)] = c1;
      en[C_LINDEX_GHOST_ONE(i,j,k)] = n1;
      hc[C_LINDEX_GHOST_ONE(i,j,k)] = c2;
      hn[C_LINDEX_GHOST_ONE(i,j,k)] = n2;
    }

  }

  restart_file.close();
}
/* -------------------------------------------------------------------------- */
/** @brief 输出粒子的浓度和数目，供重启动用.
 * */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::output_for_restart() {
  int i,j,k;
  Epetra_Vector c_ele_charge(*c_map);
  Epetra_Vector c_hole_charge(*c_map);
  Epetra_IntVector c_hole_num(*c_map);
  Epetra_IntVector c_ele_num(*c_map);
  double * ec, * hc;
  int * h_num, * e_num;
  ofstream ofile;
  int tmp;
  MPI_Status status;
  string filename;

  list<Particle> * c_par_list;

  list<Particle>::iterator iter;

  c_ele_charge.ExtractView(&ec);
  c_hole_charge.ExtractView(&hc);
  c_ele_num.ExtractView(&e_num);
  c_hole_num.ExtractView(&h_num);

  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)

    /*notice that the index starts from c_jbegin_ghost instead of c_jbegin,
     * because we want to compute the nodes with j = p_jbegin , and ghost
     * cells with j = c_jbegin_ghost have contribution on those nodes.*/

      for (j = c_jbegin; j <= c_jend; j ++){
      /*loop the cells' particle list */
	c_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];
	ec[C_LINDEX_ONE(i,j,k)] = 0;
	hc[C_LINDEX_ONE(i,j,k)] = 0;
	e_num[C_LINDEX_ONE(i,j,k)] = 0;
	h_num[C_LINDEX_ONE(i,j,k)] = 0;

	for (iter = c_par_list->begin(); iter != c_par_list->end(); iter ++){
	  charge = iter->charge;
	  par_type = iter->par_type;

	  if (par_type == PELEC){
	    ec[C_LINDEX_ONE(i,j,k)] += fabs(charge);
	    e_num[C_LINDEX_ONE(i,j,k)] ++;
	  } else {
	    hc[C_LINDEX_ONE(i,j,k)] += fabs(charge);
	    h_num[C_LINDEX_ONE(i,j,k)] ++;
	  }
	}
      }
  if (mpi_rank != 0)
    MPI_Recv(&tmp, 1, MPI_INT, mpi_rank - 1, 99 , MPI_COMM_WORLD, &status);

  filename = getFileName("./data/par_info_restart", step);
  ofile.open(filename.c_str(), iostream::app);

  for (i = c_ibegin; i <= c_iend; i ++)
    for (j = c_jbegin; j <= c_jend; j ++)
      for (k = c_kbegin; k <= c_kend; k ++)
      ofile << i << ' ' << j << ' ' << k << ' ' << ec[C_LINDEX_ONE(i,j,k)] << ' ' << e_num[C_LINDEX_ONE(i,j,k)] << ' ' 
	<< hc[C_LINDEX_ONE(i,j,k)] << ' ' << h_num[C_LINDEX_ONE(i,j,k)] << endl;

  if (mpi_rank == mpi_size - 1)
    ofile << "-1 -1 -1" << endl;

  ofile.close();

  if (mpi_rank != mpi_size - 1)
    MPI_Send(&tmp, 1, MPI_INT, mpi_rank + 1, 99, MPI_COMM_WORLD);

}
/* -------------------------------------------------------------------------- */
/** @brief 输出统计的结果.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::output_stat() {

  double * vxE, *vxH, *vzE, *vzH,*vyE, *vyH, * chargeE , * chargeH, *energyE, *energyH;
  double * stat_pot_val, * vol;
  double * stat_qc_pot_val, * save_pot_val, * save_qc_pot_val;
  double Charge, fabsCharge;
  double * saved_e_charge, * saved_h_charge;
  ofstream ofile;
  int i,j,k;
  MPI_Status status;
  string filename;
  double * c_enum, *e_heat, *h_heat;
  
  stat_pot->ExtractView(&stat_pot_val);
  stat_qc_pot->ExtractView(&stat_qc_pot_val);

  stat_pot_saved->ExtractView(&save_pot_val);
  stat_qc_pot_saved->ExtractView(&save_qc_pot_val);

  stat_vxE->ExtractView(&vxE);
  stat_vxH->ExtractView(&vxH);
  stat_vyE->ExtractView(&vyE);
  stat_vyH->ExtractView(&vyH);
  stat_vzE->ExtractView(&vzE);
  stat_vzH->ExtractView(&vzH);

  stat_chargeE->ExtractView(&chargeE);
  stat_chargeH->ExtractView(&chargeH);
  stat_energyE->ExtractView(&energyE);
  stat_energyH->ExtractView(&energyH);

  stat_chargeE_saved->ExtractView(&saved_e_charge);
  stat_chargeH_saved->ExtractView(&saved_h_charge);

  p_volume->ExtractView(&vol);

  stat_enum->ExtractView(&c_enum);

  current_scatter_info();

  if (mpi_rank != 0)
    MPI_Recv(&j, 1, MPI_INT, mpi_rank - 1, 99 , MPI_COMM_WORLD, &status);

  if (flag_heat) {
   filename = getFileName("./data/heat", step);
   ofile.open(filename.c_str(), iostream::app);

   stat_e_heat->ExtractView(&e_heat);
   stat_h_heat->ExtractView(&h_heat);

   double val1 , val2;

   for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
	{

          if (fabs(chargeE[P_LINDEX_ONE(i,j,k)]) > MY_ZERO)
            val1 = e_heat[P_LINDEX_ONE(i,j,k)] * pot0 / chargeE[P_LINDEX_ONE(i,j,k)];
          else val1 = 0;

          if (fabs(chargeH[P_LINDEX_ONE(i,j,k)]) > MY_ZERO)
	    val2 = h_heat[P_LINDEX_ONE(i,j,k)] * pot0 / chargeH[P_LINDEX_ONE(i,j,k)];
          else val2 = 0;

	  ofile << i << ' ' << j << ' ' << k << ' ' << val1  << ' ' << val1 * saved_e_charge[P_LINDEX_ONE(i,j,k)] << ' ' << val2 << ' ' << val2 * saved_h_charge[P_LINDEX_ONE(i,j,k)] << endl;
	}

   ofile.close();
  }
  
  filename = getFileName("./data/electron_num", step);
  ofile.open(filename.c_str(), iostream::app);

  for (i = c_ibegin; i <= c_iend; i ++)
    for (j = c_jbegin; j <= c_jend; j ++)
      for (k = c_kbegin; k <= c_kend; k ++)
      ofile << i << ' ' << j << ' ' << k << ' ' << c_enum[C_LINDEX_GHOST_ONE(i,j,k)] / stat_step  << endl;

  ofile.close();

  filename = getFileName("./data/pot", step);
  ofile.open(filename.c_str(), iostream::app);

  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
	{
          save_pot_val[P_LINDEX_ONE(i,j,k)] = stat_pot_val[P_LINDEX_ONE(i,j,k)] / stat_step;

	  ofile << i << ' ' << j << ' ' << k << ' ' << save_pot_val[P_LINDEX_ONE(i,j,k)] * pot0 << endl;
	}

  ofile.close();

  filename = getFileName("./data/qc_pot", step);
  ofile.open(filename.c_str(), iostream::app);

  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      {
	save_qc_pot_val[P_LINDEX_ONE(i,j,k)] = stat_qc_pot_val[P_LINDEX_ONE(i,j,k)] / stat_step;

	ofile << i << ' ' << j << ' ' << k << ' ' << save_qc_pot_val[P_LINDEX_ONE(i,j,k)] * pot0 << endl;
      }

  ofile.close();

  filename = getFileName("./data/Electron", step);
  ofile.open(filename.c_str(), iostream::app);
  
  
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
    for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++){
      Charge = chargeE[P_LINDEX_ONE(i,j,k)];

      fabsCharge = fabs(Charge);
      if (fabsCharge > MY_ZERO){
        vxE[P_LINDEX_ONE(i,j,k)] = vxE[P_LINDEX_ONE(i,j,k)] / fabsCharge * velo0;
        vyE[P_LINDEX_ONE(i,j,k)] = vyE[P_LINDEX_ONE(i,j,k)] / fabsCharge * velo0;
        vzE[P_LINDEX_ONE(i,j,k)] = vzE[P_LINDEX_ONE(i,j,k)] / fabsCharge * velo0;
        energyE[P_LINDEX_ONE(i,j,k)] = energyE[P_LINDEX_ONE(i,j,k)] / fabsCharge * pot0;
      } else{
        vxE[P_LINDEX_ONE(i,j,k)] = 0;
        vyE[P_LINDEX_ONE(i,j,k)] = 0;
        vzE[P_LINDEX_ONE(i,j,k)] = 0;
        energyE[P_LINDEX_ONE(i,j,k)] = 0;
      }
      
      if (vol[P_LINDEX_ONE(i,j,k)] > MY_ZERO)
        Charge = Charge / vol[P_LINDEX_ONE(i,j,k)] / stat_step * conc0; else Charge = 0;

      saved_e_charge[P_LINDEX_ONE(i,j,k)] = Charge;
      
      ofile << i << ' ' << j << ' ' << k << ' '
	    << vxE[P_LINDEX_ONE(i,j,k)] << ' '
            << vyE[P_LINDEX_ONE(i,j,k)] << ' '
            << vzE[P_LINDEX_ONE(i,j,k)] << ' '
            << energyE[P_LINDEX_ONE(i,j,k)] << ' '
            << Charge << endl;
    }
  
  ofile.close();
  
  if (mpi_rank != mpi_size - 1)
    MPI_Send(&j, 1, MPI_INT, mpi_rank + 1, 99, MPI_COMM_WORLD);

  if (mpi_rank != 0)
    MPI_Recv(&j, 1, MPI_INT, mpi_rank - 1, 99 , MPI_COMM_WORLD, &status);

  filename = getFileName("./data/Hole", step) ;
  ofile.open(filename.c_str(), iostream::app);
  
  
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++){
      Charge = chargeH[P_LINDEX_ONE(i,j,k)];
      fabsCharge = fabs(Charge);
      if (fabsCharge > MY_ZERO){
        vxH[P_LINDEX_ONE(i,j,k)] = vxH[P_LINDEX_ONE(i,j,k)] / fabsCharge * velo0;
        vyH[P_LINDEX_ONE(i,j,k)] = vyH[P_LINDEX_ONE(i,j,k)] / fabsCharge * velo0;
        vzH[P_LINDEX_ONE(i,j,k)] = vzH[P_LINDEX_ONE(i,j,k)] / fabsCharge * velo0;
        energyH[P_LINDEX_ONE(i,j,k)] = energyH[P_LINDEX_ONE(i,j,k)] / fabsCharge * pot0;
      } else{
        vxH[P_LINDEX_ONE(i,j,k)] = 0;
        vyH[P_LINDEX_ONE(i,j,k)] = 0;
        vzH[P_LINDEX_ONE(i,j,k)] = 0;
        energyH[P_LINDEX_ONE(i,j,k)] = 0;
      }

      if (vol[P_LINDEX_ONE(i,j,k)] > MY_ZERO)
        Charge = Charge / vol[P_LINDEX_ONE(i,j,k)] / stat_step * conc0; else Charge = 0;

      saved_h_charge[P_LINDEX_ONE(i,j,k)] = Charge;

      ofile << i << ' ' << j << ' ' << k <<' ' 
	    << vxH[P_LINDEX_ONE(i,j,k)] << ' '
            << vyH[P_LINDEX_ONE(i,j,k)] << ' '
            << vzH[P_LINDEX_ONE(i,j,k)] << ' '
            << energyH[P_LINDEX_ONE(i,j,k)] << ' '
            << Charge << endl;
    }
  
  ofile.close();
  
  if (mpi_rank != mpi_size - 1)
    MPI_Send(&j, 1, MPI_INT, mpi_rank + 1, 99, MPI_COMM_WORLD);
}
/* -------------------------------------------------------------------------- */
/** @brief 每次模拟后统计.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::statistic() {

  int i,j,k;
  double ratio_x, ratio_y, ratio_z;
  double * vxE, *vxH, *vyE, *vyH, *vzE, *vzH, * chargeE , * chargeH, *energyE, *energyH;
  double * qc_pot, * pot, *stat_pot_val, *stat_qc_pot_val, *stat_e_heat_val, * stat_h_heat_val;
  double cc, fabscc, *cvol, *e_heat, *h_heat;
  double * tmp, *ec;
  double *c_enum;

  list<Particle> * c_par_list;

  list<Particle>::iterator iter;

  stat_pot->ExtractView(&stat_pot_val);
  stat_qc_pot->ExtractView(&stat_qc_pot_val);

  /*the value get from the poisson solver, without quantum correction*/
  p_poisson_pot_saved->ExtractView(&pot);
  p_poisson_pot->ExtractView(&qc_pot);

  stat_e_heat->ExtractView(&stat_e_heat_val);
  stat_h_heat->ExtractView(&stat_h_heat_val);

  p_electron_heat->ExtractView(&e_heat);
  p_hole_heat->ExtractView(&h_heat);

  /*get the summation for each node */
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++){
	      stat_pot_val[P_LINDEX_ONE(i,j,k)] += pot[P_LINDEX_ONE(i,j,k)];
	      stat_qc_pot_val[P_LINDEX_ONE(i,j,k)] += qc_pot[P_LINDEX_ONE(i,j,k)];
        stat_e_heat_val[P_LINDEX_ONE(i,j,k)] += e_heat[P_LINDEX_ONE_GHOST(i,j,k)];
        stat_h_heat_val[P_LINDEX_ONE(i,j,k)] += h_heat[P_LINDEX_ONE_GHOST(i,j,k)];
    }

  /*prepare data array*/ 
  stat_vxE->ExtractView(&vxE);
  stat_vxH->ExtractView(&vxH);
  stat_vyE->ExtractView(&vyE);
  stat_vyH->ExtractView(&vyH);
  stat_vzE->ExtractView(&vzE);
  stat_vzH->ExtractView(&vzH);

  stat_chargeE->ExtractView(&chargeE);
  stat_chargeH->ExtractView(&chargeH);
  stat_energyE->ExtractView(&energyE);
  stat_energyH->ExtractView(&energyH);
  stat_enum->ExtractView(&c_enum);

  stat_ec->ExtractView(&ec);

  c_volume->ExtractView(&cvol);

  /*loop for each cell */
  for (i = c_ibegin; i <= c_iend; i ++)

    for (k = c_kbegin; k <= c_kend; k ++)

    /*notice that the index starts from c_jbegin_ghost instead of c_jbegin,
     * because we want to compute the nodes with j = p_jbegin , and ghost
     * cells with j = c_jbegin_ghost have contribution on those nodes.*/

      for (j = c_jbegin_ghost; j <= c_jend; j ++){
      
      /*loop the cells' particle list */
	c_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];
        c_enum[C_LINDEX_GHOST_ONE(i,j,k)] += c_par_list->size();

//      if (NOT_GHOST_CELL(j)){
      {
      for (iter = c_par_list->begin(); iter != c_par_list->end(); iter ++){
	/*fill in the static variable of the class */
        InPar(iter);
	/* get velocity based on k-state */
        GetV();
	/* the particle's contribution to each node is based on the position
	 * in the cell, roughly speaking, the further away the particle is from the node, the 
	 * less the contribution is.*/
        ratio_x = (x - lx[i]) / dx[i];
        ratio_y = (y - ly[j]) / dy[j];
        ratio_z = (z - lz[k]) / dz[k];

	/*simple assertion: particle should be in this cell*/
	if (!(BETWEEN01(ratio_x) && BETWEEN01(ratio_y) && BETWEEN01(ratio_z))) {
	  err_message(WRONG_CELL, "statistic");
	  dump_par_info();
	}

        /*for electrons*/ 
        if (par_type == PELEC){

	  ec[C_LINDEX_GHOST_ONE(i,j,k)] += fabs(charge) / cvol[C_LINDEX_GHOST_ONE(i,j,k)] / stat_step ;
          

          if (j >= p_jbegin_nonoverlap){  /*contribution on the left side node*/
            cc = charge * (1 - ratio_x) * (1 - ratio_y) * (1 - ratio_z);
            fabscc = fabs(cc);
            vxE[P_LINDEX_ONE(i,j,k)] += vx * fabscc;
            vyE[P_LINDEX_ONE(i,j,k)] += vy * fabscc;
            vzE[P_LINDEX_ONE(i,j,k)] += vz * fabscc;
            chargeE[P_LINDEX_ONE(i,j,k)] += fabscc;
            energyE[P_LINDEX_ONE(i,j,k)] += fabscc * energy;

            cc = charge * ratio_x * (1 - ratio_y) * (1 - ratio_z);
            fabscc = fabs(cc);
            vxE[P_LINDEX_ONE(i + 1,j,k)] += vx * fabscc;
            vyE[P_LINDEX_ONE(i + 1,j,k)] += vy * fabscc;
            vzE[P_LINDEX_ONE(i + 1,j,k)] += vz * fabscc;
            chargeE[P_LINDEX_ONE(i + 1,j,k)] += fabscc;
            energyE[P_LINDEX_ONE(i + 1,j,k)] += fabscc * energy;

	    cc = charge * (1 - ratio_x) * (1 - ratio_y) * ratio_z;
            fabscc = fabs(cc);
            vxE[P_LINDEX_ONE(i,j,k + 1)] += vx * fabscc;
            vyE[P_LINDEX_ONE(i,j,k+1)] += vy * fabscc;
            vzE[P_LINDEX_ONE(i,j,k+1)] += vz * fabscc;
            chargeE[P_LINDEX_ONE(i,j,k+1)] += fabscc;
            energyE[P_LINDEX_ONE(i,j,k+1)] += fabscc * energy;

            cc = charge * ratio_x * (1 - ratio_y) * ratio_z;
            fabscc = fabs(cc);
            vxE[P_LINDEX_ONE(i + 1, j, k+1)] += vx * fabscc;
            vyE[P_LINDEX_ONE(i + 1, j, k+1)] += vy * fabscc;
            vzE[P_LINDEX_ONE(i + 1, j, k+1)] += vz * fabscc;
            chargeE[P_LINDEX_ONE(i + 1,j,k+1)] += fabscc;
            energyE[P_LINDEX_ONE(i + 1,j,k+1)] += fabscc * energy;

          }
	  /*notice that when j = c_jend, we need to check 
	   * whether the right side node should be computed to avoid the index
	   * go out of the range */
	  if (j < p_jend_nonoverlap) {    /* contribution on the right side node */
	    cc = charge * (1 - ratio_x) * ratio_y * (1 - ratio_z);
	    fabscc = fabs(cc);
	    vxE[P_LINDEX_ONE(i,j + 1, k)] += vx * fabscc;
	    vyE[P_LINDEX_ONE(i,j + 1,k)] += vy * fabscc;
	    vzE[P_LINDEX_ONE(i,j + 1,k)] += vz * fabscc;
	    chargeE[P_LINDEX_ONE(i,j + 1,k)] += fabscc;
	    energyE[P_LINDEX_ONE(i,j + 1,k)] += fabscc * energy;
	    
	    cc = charge * ratio_x * ratio_y * (1 - ratio_z);
	    fabscc = fabs(cc);
	    vxE[P_LINDEX_ONE(i + 1,j + 1,k)] += vx * fabscc;
	    vyE[P_LINDEX_ONE(i + 1,j + 1,k)] += vy * fabscc;
	    vzE[P_LINDEX_ONE(i + 1,j + 1,k)] += vz * fabscc;
	    chargeE[P_LINDEX_ONE(i + 1,j + 1,k)] += fabscc;
	    energyE[P_LINDEX_ONE(i + 1,j + 1,k)] += fabscc * energy;

	    cc = charge * (1 - ratio_x) * ratio_y * ratio_z;
	    fabscc = fabs(cc);
	    vxE[P_LINDEX_ONE(i,j + 1, k+1)] += vx * fabscc;
	    vyE[P_LINDEX_ONE(i,j + 1,k+1)] += vy * fabscc;
	    vzE[P_LINDEX_ONE(i,j + 1,k+1)] += vz * fabscc;
	    chargeE[P_LINDEX_ONE(i,j + 1,k+1)] += fabscc;
	    energyE[P_LINDEX_ONE(i,j + 1,k+1)] += fabscc * energy;
	    
	    cc = charge * ratio_x * ratio_y * ratio_z;
	    fabscc = fabs(cc);
	    vxE[P_LINDEX_ONE(i + 1,j + 1,k+1)] += vx * fabscc;
	    vyE[P_LINDEX_ONE(i + 1,j + 1,k+1)] += vy * fabscc;
	    vzE[P_LINDEX_ONE(i + 1,j + 1,k+1)] += vz * fabscc;
	    chargeE[P_LINDEX_ONE(i + 1,j + 1,k+1)] += fabscc;
	    energyE[P_LINDEX_ONE(i + 1,j + 1,k+1)] += fabscc * energy;
	  }
	}

        /*for holes */
        if (par_type == PHOLE){
          if (j >= p_jbegin_nonoverlap){/*contribution on the left side node*/
            cc = charge * (1 - ratio_x) * (1 - ratio_y) * (1 - ratio_z);
            fabscc = fabs(cc);
            vxH[P_LINDEX_ONE(i,j,k)] += vx * fabscc;
            vyH[P_LINDEX_ONE(i,j,k)] += vy * fabscc;
            vzH[P_LINDEX_ONE(i,j,k)] += vz * fabscc;
            chargeH[P_LINDEX_ONE(i,j,k)] += fabscc;
            energyH[P_LINDEX_ONE(i,j,k)] += fabscc * energy;

            cc = charge * ratio_x * (1 - ratio_y)* (1 - ratio_z);
            fabscc = fabs(cc);
            vxH[P_LINDEX_ONE(i + 1,j,k)] += vx * fabscc;
            vyH[P_LINDEX_ONE(i + 1,j,k)] += vy * fabscc;
            vzH[P_LINDEX_ONE(i + 1,j,k)] += vz * fabscc;
            chargeH[P_LINDEX_ONE(i + 1,j,k)] += fabscc;
            energyH[P_LINDEX_ONE(i + 1,j,k)] += fabscc * energy;

            cc = charge * (1 - ratio_x) * (1 - ratio_y) * ratio_z;
            fabscc = fabs(cc);
            vxH[P_LINDEX_ONE(i,j,k+1)] += vx * fabscc;
            vyH[P_LINDEX_ONE(i,j,k+1)] += vy * fabscc;
            vzH[P_LINDEX_ONE(i,j,k+1)] += vz * fabscc;
            chargeH[P_LINDEX_ONE(i,j,k+1)] += fabscc;
            energyH[P_LINDEX_ONE(i,j,k+1)] += fabscc * energy;

            cc = charge * ratio_x * (1 - ratio_y) * ratio_z;
            fabscc = fabs(cc);
            vxH[P_LINDEX_ONE(i + 1,j,k+1)] += vx * fabscc;
            vyH[P_LINDEX_ONE(i + 1,j,k+1)] += vy * fabscc;
            vzH[P_LINDEX_ONE(i + 1,j,k+1)] += vz * fabscc;
            chargeH[P_LINDEX_ONE(i + 1,j,k+1)] += fabscc;
            energyH[P_LINDEX_ONE(i + 1,j,k+1)] += fabscc * energy;

          }
          if (j < p_jend_nonoverlap) { /*contribution on the right side node */
	    cc = charge * (1 - ratio_x) * ratio_y * (1 - ratio_z);
	    fabscc = fabs(cc);
	    vxH[P_LINDEX_ONE(i,j + 1,k)] += vx * fabscc;
	    vyH[P_LINDEX_ONE(i,j + 1,k)] += vy * fabscc;
	    vzH[P_LINDEX_ONE(i,j + 1,k)] += vz * fabscc;
	    chargeH[P_LINDEX_ONE(i,j + 1,k)] += fabscc;
	    energyH[P_LINDEX_ONE(i,j + 1,k)] += fabscc * energy;
	    
	    cc = charge * ratio_x * ratio_y * (1 - ratio_z);
	    fabscc = fabs(cc);
	    vxH[P_LINDEX_ONE(i + 1,j + 1,k)] += vx * fabscc;
	    vyH[P_LINDEX_ONE(i + 1,j + 1,k)] += vy * fabscc;
	    vzH[P_LINDEX_ONE(i + 1,j + 1,k)] += vz * fabscc;
	    chargeH[P_LINDEX_ONE(i + 1,j + 1,k)] += fabscc;
	    energyH[P_LINDEX_ONE(i + 1,j + 1,k)] += fabscc * energy;

	    cc = charge * (1 - ratio_x) * ratio_y * ratio_z;
	    fabscc = fabs(cc);
	    vxH[P_LINDEX_ONE(i,j + 1,k+1)] += vx * fabscc;
	    vyH[P_LINDEX_ONE(i,j + 1,k+1)] += vy * fabscc;
	    vzH[P_LINDEX_ONE(i,j + 1,k+1)] += vz * fabscc;
	    chargeH[P_LINDEX_ONE(i,j + 1,k+1)] += fabscc;
	    energyH[P_LINDEX_ONE(i,j + 1,k+1)] += fabscc * energy;
	    
	    cc = charge * ratio_x * ratio_y * ratio_z;
	    fabscc = fabs(cc);
	    vxH[P_LINDEX_ONE(i + 1,j + 1,k+1)] += vx * fabscc;
	    vyH[P_LINDEX_ONE(i + 1,j + 1,k+1)] += vy * fabscc;
	    vzH[P_LINDEX_ONE(i + 1,j + 1,k+1)] += vz * fabscc;
	    chargeH[P_LINDEX_ONE(i + 1,j + 1,k+1)] += fabscc;
	    energyH[P_LINDEX_ONE(i + 1,j + 1,k+1)] += fabscc * energy;
	  }
	}
      }
      }
    }
}
/* -------------------------------------------------------------------------- */
/** @brief 检查粒子数是否守恒，防止因程序 bug 丢粒子
 * */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::check_par_number() {

  int tot_catch_par, tot_gen_par;
  int local_actual_par, global_actual_par, tot_mr_gen_num;
  int i,j,k;
  
  MPI_Allreduce(&catch_par_num, &tot_catch_par, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&gen_par, &tot_gen_par, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mr_gen_num, &tot_mr_gen_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  local_actual_par = 0;

  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
      for (j = c_jbegin; j <= c_jend; j ++){
	local_actual_par += par_list[C_LINDEX_GHOST_ONE(i,j,k)].size();
      }

  MPI_Allreduce(&local_actual_par, &global_actual_par, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  if (mpi_rank == 0){

    if (global_actual_par != par_num - tot_catch_par + tot_gen_par + tot_mr_gen_num){
       cout << "particle number wrong: " << endl
//             << "inject particles " << inject_par_num<< endl
             << "catch " << tot_catch_par << endl
             << "generate " << tot_gen_par<< endl
	     << "mr_gen " << tot_mr_gen_num << endl
             << "expected " << par_num - tot_catch_par + tot_gen_par + tot_mr_gen_num << endl
             << "actual " << global_actual_par << endl;
      exit(1);
    }
  }
  
  par_num = global_actual_par;
  
  if (mpi_rank == 0) {
    cout << " generate : " << tot_gen_par  
         << ", catch : " << tot_catch_par << endl;
  }
}
/* -------------------------------------------------------------------------- */
/** @brief 模拟粒子的飞行过程
 */
/* ---------------------------------------------------------------------------- */
void MeshQuantities::update_particle() {
  int i,j,k;
  int unfinish_par_tot = 0, unfinish_par;
  int loop = 0;

  catch_par_num = 0;
  gen_par = 0;

  if (mpi_rank == 0)
    cout << "particle fly .. " << endl;

  do {
    
    loop ++;
    
    unfinish_par = 0;
    
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	      for (j = c_jbegin; j <= c_jend; j ++){
	        current_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];

          /*	if ((!NOT_GHOST_CELL(j)) && (current_par_list->size() != 0))
	          err_message(GHOST_PARTICLE, "should not have particles in ghosts cells");
	        */

          for (par_iter = current_par_list->begin(); par_iter != current_par_list->end(); par_iter ++)
            if ((loop == 1) || (par_iter->flag == 0)) {

              particle_fly();
            
              if ((par_iter->i >= 0) && (par_iter->flag == 0)) 
                unfinish_par ++;
            }
        }
    
    MPI_Allreduce(&unfinish_par, &unfinish_par_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (mpi_rank == 0)
      cout << " loop " << loop << " : " << unfinish_par_tot << " particles unfinished " << endl;
    
    migrate_particle();
    
  } while (unfinish_par_tot > 0);

}
/* -------------------------------------------------------------------------- */
/** @brief 把粒子迁移到正确的位置.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::migrate_particle() {
  
  local_migrate();
  
  global_migrate();

  clear_ghost_par();
}

void MeshQuantities::fill_ghost_cell() {

  int par_send_num, par_recv_num;
  MPI_Request request, request_int;
  MPI_Status status;
  
  double *dbl_snd_buf, *dbl_snd_buf_right, *dbl_recv_buf, *dbl_recv_buf_right;
  int *int_snd_buf, *int_snd_buf_right, *int_recv_buf, *int_recv_buf_right;
  int i,k;
  par_send_num = 0;
  par_recv_num = 0;

  if (mpi_rank != mpi_size - 1) {
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	par_send_num += par_list[C_LINDEX_GHOST_ONE(i, c_jend,k)].size();
    
    MPI_Isend(&par_send_num, 1, MPI_INT, mpi_rank + 1, 1, MPI_COMM_WORLD, &request);
  }

  if (mpi_rank != 0) 
    MPI_Recv(&par_recv_num, 1, MPI_INT, mpi_rank - 1, 1, MPI_COMM_WORLD, &status);
  
  if (mpi_rank != mpi_size - 1)
    MPI_Wait(&request, &status);

  /*if we will send some particles to mpi_rank + 1 */
  if ((mpi_rank != mpi_size - 1) && (par_send_num > 0)) {
    int l = 0;
    dbl_snd_buf =(double *) malloc(sizeof(double) * PAR_DBL_NUM * par_send_num);
    int_snd_buf = (int *) malloc(sizeof(int) * PAR_INT_NUM * par_send_num);
    
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	for (par_iter = par_list[C_LINDEX_GHOST_ONE(i,c_jend,k)].begin(); 
           par_iter != par_list[C_LINDEX_GHOST_ONE(i,c_jend,k)].end(); par_iter ++){
	  cp_to_buf(&dbl_snd_buf[PAR_DBL_NUM * l],&int_snd_buf[PAR_INT_NUM * l], *par_iter);
	  l ++;
	}
    MPI_Isend(dbl_snd_buf, PAR_DBL_NUM * par_send_num, MPI_DOUBLE, mpi_rank + 1, 1, MPI_COMM_WORLD, &request);
    MPI_Isend(int_snd_buf, PAR_INT_NUM * par_send_num, MPI_INT, mpi_rank + 1, 2, MPI_COMM_WORLD, &request_int);
  }
  /*if we will receive some particles from mpi_rank - 1 */
  if ((mpi_rank != 0) && (par_recv_num > 0)) {
    dbl_recv_buf =(double *) malloc(sizeof(double) * PAR_DBL_NUM * par_recv_num);
    int_recv_buf = (int *) malloc(sizeof(int) * PAR_INT_NUM * par_recv_num);
    MPI_Recv(dbl_recv_buf, PAR_DBL_NUM * par_recv_num, MPI_DOUBLE, mpi_rank - 1, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(int_recv_buf, PAR_INT_NUM * par_recv_num, MPI_INT, mpi_rank - 1, 2, MPI_COMM_WORLD, &status);
  }
  if ((mpi_rank != mpi_size - 1) && (par_send_num > 0)){
    MPI_Wait(&request, &status);
    MPI_Wait(&request_int, &status);
  }
  
  Particle par_tmp;
  
  if (mpi_rank != 0) {
    for (i = 0;i < par_recv_num; i ++){
      par_tmp = make_particle(&dbl_recv_buf[PAR_DBL_NUM * i], &int_recv_buf[PAR_INT_NUM * i]);
      par_list[C_LINDEX_GHOST_ONE(par_tmp.i, par_tmp.j, par_tmp.k)].push_back(par_tmp);
    }
  }
  
  if ((mpi_rank != mpi_size - 1) && (par_send_num > 0)) {
    free(dbl_snd_buf);
    free(int_snd_buf);
  }
  if ((mpi_rank != 0) && (par_recv_num > 0)) {
    free(dbl_recv_buf);
    free(int_recv_buf);
  }
}
/* -------------------------------------------------------------------------- */
/** @brief 清楚 ghost 单元格中的粒子.
 */
/* ---------------------------------------------------------------------------- */
 
void MeshQuantities::clear_ghost_par() {
  int i,k;
  if (mpi_rank != 0) {
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	par_list[C_LINDEX_GHOST_ONE(i, c_jbegin_ghost,k)].clear();
  }
  if (mpi_rank != mpi_size - 1) {
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	par_list[C_LINDEX_GHOST_ONE(i, c_jend_ghost,k)].clear();
  }
}
/* -------------------------------------------------------------------------- */
/** @brief 在当前节点内做迁移.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::local_migrate() {
  int i,j,k;
  
  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
    
      for (j = c_jbegin; j <= c_jend; j ++){
	
	current_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];
	
	for (par_iter = current_par_list->begin(); par_iter != current_par_list->end(); ) {

          // 先拦截负索引（正常 catch 会设为 -9999，其他负值视为异常逃逸）
          if (par_iter->i < 0) {
              if (par_iter->i != -9999) {
                  cout << "[CRITICAL LEAK] Particle deleted with negative index!" << endl;
                  cout << "  ID: " << par_iter->par_id << endl;
                  cout << "  Index i: " << par_iter->i << " (Should be -9999 if caught)" << endl;
                  cout << "  Pos: " << par_iter->x * spr0 << ", " << par_iter->y * spr0 << ", " << par_iter->z * spr0 << endl;
                  cout << "  From Cell: (" << i << "," << j << "," << k << ")" << endl;
              }
              par_iter = current_par_list->erase(par_iter);
              continue;
          }

	  if ((par_iter->i != i) || (par_iter->j != j) || (par_iter->k != k)) {

            // Debug: 粒子是否越界 (X/Y/Z)
            bool is_j_out = (par_iter->j < c_jbegin_ghost || par_iter->j > c_jend_ghost);
            bool is_i_out = (par_iter->i < c_ibegin || par_iter->i > c_iend);
            bool is_k_out = (par_iter->k < c_kbegin || par_iter->k > c_kend);

            if (is_j_out || is_i_out || is_k_out) {
                cout << "[CRITICAL LEAK] Particle Escaped Domain!" << endl;
                cout << "  ID: " << par_iter->par_id << " Type: " << par_iter->par_type << endl;
                cout << "  From Cell: (" << i << ", " << j << ", " << k << ")" << endl;
                cout << "  To Cell  : (" << par_iter->i << ", " << par_iter->j << ", " << par_iter->k << ")" << endl;
                cout << "  Position : (" << par_iter->x * spr0 << ", " << par_iter->y * spr0 << ", " << par_iter->z * spr0 << ")" << endl;
                
                if (is_i_out) cout << "  -> X-direction LEAK (i range: " << c_ibegin << "~" << c_iend << ")" << endl;
                if (is_j_out) cout << "  -> Y-direction LEAK (j range: " << c_jbegin_ghost << "~" << c_jend_ghost << ")" << endl;
                if (is_k_out) cout << "  -> Z-direction LEAK (k range: " << c_kbegin << "~" << c_kend << ")" << endl;

                // 越界粒子直接删除，避免 push_back 越界
                par_iter = current_par_list->erase(par_iter);
                continue;
            }

	    par_list[C_LINDEX_GHOST_ONE(par_iter->i, par_iter->j, par_iter->k)].push_back(*par_iter);
	    
	    par_iter = current_par_list->erase(par_iter);
	  }
	  else par_iter ++;
        }
      
    }
}
/* -------------------------------------------------------------------------- */
/** @brief 跨节点迁移.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::global_migrate() {
  int par_send_num[2], par_recv_num[2];
  MPI_Request request[2], request_int[2];
  MPI_Status status;
  
  double *dbl_snd_buf, *dbl_snd_buf_right, *dbl_recv_buf, *dbl_recv_buf_right;
  int *int_snd_buf, *int_snd_buf_right, *int_recv_buf, *int_recv_buf_right;
  int i,k;
  par_send_num[0] = 0;
  par_send_num[1] = 0;
  par_recv_num[0] = 0;
  par_recv_num[1] = 0;

  if (mpi_rank != 0) {
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	par_send_num[0] += par_list[C_LINDEX_GHOST_ONE(i, c_jbegin_ghost, k)].size();
    MPI_Isend(&par_send_num[0], 1, MPI_INT, mpi_rank - 1, 1, MPI_COMM_WORLD, &request[0]);
  }
  if (mpi_rank != mpi_size - 1) {
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
      par_send_num[1] += par_list[C_LINDEX_GHOST_ONE(i, c_jend_ghost,k)].size();
    
    MPI_Isend(&par_send_num[1], 1, MPI_INT, mpi_rank + 1, 1, MPI_COMM_WORLD, &request[1]);
  }

  if (mpi_rank != 0) 
    MPI_Recv(&par_recv_num[0], 1, MPI_INT, mpi_rank - 1, 1, MPI_COMM_WORLD, &status);
  
  if (mpi_rank != mpi_size - 1) 
    MPI_Recv(&par_recv_num[1], 1, MPI_INT, mpi_rank + 1, 1, MPI_COMM_WORLD, &status);
  
  if (mpi_rank != 0)
    MPI_Wait(&request[0], &status);
  
  if (mpi_rank != mpi_size - 1)
    MPI_Wait(&request[1], &status);

  if ((mpi_rank != 0) && (par_send_num[0] > 0)) {
    int l = 0;
    dbl_snd_buf =(double *) malloc(sizeof(double) * PAR_DBL_NUM * par_send_num[0]);
    int_snd_buf = (int *) malloc(sizeof(int) * PAR_INT_NUM * par_send_num[0]);
      
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
      for (par_iter = par_list[C_LINDEX_GHOST_ONE(i,c_jbegin_ghost,k)].begin(); 
             par_iter != par_list[C_LINDEX_GHOST_ONE(i,c_jbegin_ghost,k)].end(); par_iter ++){
        cp_to_buf(&dbl_snd_buf[PAR_DBL_NUM * l], &int_snd_buf[PAR_INT_NUM * l],  *par_iter);
        l ++;
      }
    MPI_Isend(dbl_snd_buf, PAR_DBL_NUM * par_send_num[0], MPI_DOUBLE, mpi_rank - 1, 1, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(int_snd_buf, PAR_INT_NUM * par_send_num[0], MPI_INT, mpi_rank - 1, 2, MPI_COMM_WORLD, &request_int[0]);
  }
  
  if ((mpi_rank != mpi_size - 1) && (par_send_num[1] > 0)) {
    int l = 0;
    dbl_snd_buf_right =(double *) malloc(sizeof(double) * PAR_DBL_NUM * par_send_num[1]);
    int_snd_buf_right = (int *) malloc(sizeof(int) * PAR_INT_NUM * par_send_num[1]);
    
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
	for (par_iter = par_list[C_LINDEX_GHOST_ONE(i,c_jend_ghost,k)].begin(); 
             par_iter != par_list[C_LINDEX_GHOST_ONE(i,c_jend_ghost,k)].end(); par_iter ++){
        cp_to_buf(&dbl_snd_buf_right[PAR_DBL_NUM * l],&int_snd_buf_right[PAR_INT_NUM * l], *par_iter);
        l ++;
      }
    MPI_Isend(dbl_snd_buf_right, PAR_DBL_NUM * par_send_num[1], MPI_DOUBLE, mpi_rank + 1, 1, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(int_snd_buf_right, PAR_INT_NUM * par_send_num[1], MPI_INT, mpi_rank + 1, 2, MPI_COMM_WORLD, &request_int[1]);
  }
  if ((mpi_rank != 0) && (par_recv_num[0] > 0)) {
    dbl_recv_buf =(double *) malloc(sizeof(double) * PAR_DBL_NUM * par_recv_num[0]);
    int_recv_buf = (int *) malloc(sizeof(int) * PAR_INT_NUM * par_recv_num[0]);
    MPI_Recv(dbl_recv_buf, PAR_DBL_NUM * par_recv_num[0], MPI_DOUBLE, mpi_rank - 1, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(int_recv_buf, PAR_INT_NUM * par_recv_num[0], MPI_INT, mpi_rank - 1, 2, MPI_COMM_WORLD, &status);
  }
  if ((mpi_rank != mpi_size - 1) && (par_recv_num[1] > 0)) {
    dbl_recv_buf_right =(double *) malloc(sizeof(double) * PAR_DBL_NUM * par_recv_num[1]);
    int_recv_buf_right = (int *) malloc(sizeof(int) * PAR_INT_NUM * par_recv_num[1]);
    MPI_Recv(dbl_recv_buf_right, PAR_DBL_NUM * par_recv_num[1], MPI_DOUBLE, mpi_rank + 1, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(int_recv_buf_right, PAR_INT_NUM * par_recv_num[1], MPI_INT, mpi_rank + 1, 2, MPI_COMM_WORLD, &status);
  }

  if ((mpi_rank != 0) && (par_send_num[0] > 0)) {
    MPI_Wait(&request[0], &status);
    MPI_Wait(&request_int[0], &status);
  }
  
  if ((mpi_rank != mpi_size - 1) && (par_send_num[1] > 0)){
    MPI_Wait(&request[1], &status);
    MPI_Wait(&request_int[1], &status);
  }
  
  Particle par_tmp;
  
  if (mpi_rank != 0) {
    for (i = 0;i < par_recv_num[0]; i ++){
      par_tmp = make_particle(&dbl_recv_buf[PAR_DBL_NUM * i], &int_recv_buf[PAR_INT_NUM * i]);
      par_list[C_LINDEX_GHOST_ONE(par_tmp.i, par_tmp.j,par_tmp.k)].push_back(par_tmp);
    }
  }
  
  if (mpi_rank != mpi_size - 1){
    for (i = 0;i < par_recv_num[1]; i ++){
      par_tmp = make_particle(&dbl_recv_buf_right[PAR_DBL_NUM * i], &int_recv_buf_right[PAR_INT_NUM * i]);
      par_list[C_LINDEX_GHOST_ONE(par_tmp.i, par_tmp.j,par_tmp.k)].push_back(par_tmp);
    }
  }
  if ((mpi_rank != 0) && (par_send_num[0] > 0)) {
    free(dbl_snd_buf);
    free(int_snd_buf);
  }
  if ((mpi_rank != mpi_size - 1) && (par_send_num[1] > 0)) {
    free(dbl_snd_buf_right);
    free(int_snd_buf_right);
  }
  if ((mpi_rank != 0) && (par_recv_num[0] > 0)) {
    free(dbl_recv_buf);
    free(int_recv_buf);
  }
  if ((mpi_rank != mpi_size - 1) && (par_recv_num[1] > 0)) {
    free(dbl_recv_buf_right);
    free(int_recv_buf_right);
  }
}
/* -------------------------------------------------------------------------- */
/** @brief 时间倒流
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::trace_back() {
    kx -= Ex * Tf;
    ky -= Ey * Tf;
    kz -= Ez * Tf;
    energy -= (vx * Ex + vy * Ey + vz * Ez) * Tf;
    x -= vx * Tf;
    y -= vy * Tf;
    z -= vz * Tf;
    left_time += Tf;
}
/* -------------------------------------------------------------------------- */
/** @brief dump 出网格的信息，调试用.
 * 
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::dump_cell_info(int i,int j, int k) {
  cout << "**********dump cell**********" << endl
       << "cell index     : " << '(' << i << ',' << j << ',' << k << ')' << endl
       << "x interval     : " << lx[i] << ',' << lx[i+ 1] << endl
       << "y interval     : " << ly[j] << ',' << ly[j+ 1] << endl
       << "z interval     : " << lz[k] << ',' << lz[k+ 1] << endl
       << "electric field : " << Ex << ' ' << Ey << ' ' << Ez << endl
       << "cell charge    : " << Rho << endl
       << "DA             : " << DA << endl
       << endl;
}
/* -------------------------------------------------------------------------- */
/** @brief dump 出时间信息，调试用.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::dump_time_info() {

   cout  <<"***********dump time**********"<<endl
	<<"tetra time          : " << TetTf << endl
	<<"cell time           : " << CellTf  << endl
	<<"phonon scatter time : " << PhScTf << endl
	<<"impact scatter time : " << ImpScTf << endl
	<<"surface scatter time: " << SurfScTf<<endl
	<< endl;
}
/* -------------------------------------------------------------------------- */
/** @brief dump 出粒子信息，调试用.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::dump_par_info() {
  cout << " *********dump particle********" << endl
       << " particle id : " << par_id << endl
       << " seed        : " << seed << endl
       << " cell index  : " << icell << ' ' << jcell << ' ' << kcell << endl
       << " x position  : " << x * spr0 <<  endl
       << " x ratio     : " << (x - lx[icell]) / dx[icell] << endl
       << " y position  : " << y * spr0 << endl
       << " y ratio     : " << (y - ly[jcell]) / dy[jcell] << endl
       << " z position  : " << z * spr0 << endl
       << " z ratio     : " << (z - lz[kcell]) / dz[kcell] << endl
       << " k vector    : " << kx << ',' << ky << ',' << kz << endl
       << " left time   : " << left_time << endl
       << " charge      : " << charge << endl
       << " energy      : " << energy << endl
       << " itet        : " << itet << endl
       << " iband       : " << iband << endl
       << " isym        : " << isym << endl
       << " velocity    : " << vx << ' ' << vy << ' ' << vz << endl
       << endl;
}
/* -------------------------------------------------------------------------- */
/** @brief dump 粒子信息，调试用. 
 * 
 * @param par 
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::dump_par_info(Particle & par) {
  cout << " *********dump particle********" << endl
       << " particle id : " << par.par_id << endl
       << " seed        : " << par.seed << endl
       << " cell index  : " << par.i<< ' ' << par.j<< ' ' << par.k<< endl
       << " x position  : " << par.x * spr0 <<  endl
       << " x ratio     : " << (par.x - lx[par.i]) / dx[par.i] << endl
       << " y position  : " << par.y  * spr0 << endl
       << " y ratio     : " << (par.y - ly[par.j]) / dy[par.j] << endl
       << " z position  : " << par.z  * spr0 << endl
       << " z ratio     : " << (par.z - lz[par.k]) / dz[par.k] << endl
       << " k vector    : " << par.kx << ',' << par.ky << ',' << par.kz << endl
       << " left time   : " << par.left_time << endl
       << " charge      : " << par.charge << endl
       << " energy      : " << par.energy << endl
       << " itet        : " << par.itet << endl
       << " iband       : " << band.ibt[par.itet]<< endl
       << " isym        : " << par.isym << endl
       << endl;
}
/* -------------------------------------------------------------------------- */
/** @brief 模拟一个粒子.
 */
/* ---------------------------------------------------------------------------- */

void MeshQuantities::particle_fly() {

  Flag_Catch = false;
	
  bool Flag_GetTetTime    = true;
  bool Flag_GetCellTime   = true;
  bool Flag_GetPhScTime   = true;
  bool Flag_GetImpScTime  = true;
  bool Flag_GetSurfScTime = true;

  double ssnl = 0, imprnl = 0, phrnl = 0;

  bool fly_too_far = false;
  
  int flag = 0, old_flag = 0;

  InPar(par_iter);
  
  Rho = (*c_par_charge)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
  DA = (*c_da)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
  
  if (par_type == 0) {
    Ex  = ChargeSign[par_type] * (*c_field_x)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    Ey  = ChargeSign[par_type] * (*c_field_y)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    Ez  = ChargeSign[par_type] * (*c_field_z)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
  } else {
    Ex  = ChargeSign[par_type] * (*c_h_field_x)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    Ey  = ChargeSign[par_type] * (*c_h_field_y)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    Ez  = ChargeSign[par_type] * (*c_h_field_z)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
  }

  iband = band.ibt[itet];

  if (band.use_analytic_band) {
      band.GetAnalyticV_FromTable(&(*par_iter));
      vx = band.analytic_vx;
      vy = band.analytic_vy;
      vz = band.analytic_vz;
  } else {
      GetV();
  }

  int loop = 0;
  bool old_flag_getTetTime;
  Particle old_par, new_par;
  int old_cellflag=0;

  // [DEBUG] 黑匣子：记录上一轮的完整状态
  struct StateSnapshot {
      int loop_idx;
      int trigger_flag; // 触发本轮的原因 (1=HitTet, 2=HitCell...)
      double x, y, z;
      double kx, ky, kz;
      double vx, vy, vz;
      int i, j, k;
      double calc_TetTf; // 上一轮计算出的 TetTf
      double calc_CellTf; // 上一轮计算出的 CellTf
  } prev_step, curr_step;

  // 初始化
  prev_step.loop_idx = 0; prev_step.trigger_flag = -1;
  prev_step.x = x; prev_step.y = y; prev_step.z = z;
  prev_step.calc_TetTf = -1; prev_step.calc_CellTf = -1;

  if(CellTime()<MY_ZERO){
 //     cout<<"debug point"<<endl;
  }

  while (left_time > 0) {
    
    // ----------------------------------------------------
    // 解析能带模式下的速度同步
    // ----------------------------------------------------
    if (band.use_analytic_band) {
        par_iter->kx = kx; par_iter->ky = ky; par_iter->kz = kz;
        par_iter->energy = energy;
        par_iter->x = x; par_iter->y = y; par_iter->z = z;
        //band.GetAnalyticV_FromTable(&(*par_iter));
        //vx = band.analytic_vx; vy = band.analytic_vy; vz = band.analytic_vz;
        //energy = par_iter->energy;
    } else {
        GetV();
    }

    loop ++;
    
    // 记录本轮开始时的状态
    curr_step.loop_idx = loop;
    curr_step.x = x; curr_step.y = y; curr_step.z = z;
    curr_step.kx = kx; curr_step.ky = ky; curr_step.kz = kz;
    curr_step.vx = vx; curr_step.vy = vy; curr_step.vz = vz;
    curr_step.i = icell; curr_step.j = jcell; curr_step.k = kcell;
    
    // --- 1. 计算时间步长 ---
    if(Flag_GetTetTime) {
        if (band.use_analytic_band) {
            TetTf = band.GetAnalyticGridTime(&(*par_iter), Ex, Ey, Ez);
        } else {
            TetTf = TetTime();    
        }
        Flag_GetTetTime = false;
    }
    
    if(Flag_GetCellTime) {
        CellTf = CellTime();   
        Flag_GetCellTime = false;
    }
    
    if(Flag_GetPhScTime) {
        phrnl = -log(Random());
        Flag_GetPhScTime = false;
    }
    PhScTf = phrnl / band.gamtet[itet];
    
    if(band.use_analytic_band || iband == band.bandof[PELEC]) {
        if(Flag_GetImpScTime) {
              imprnl = -log(Random());
              Flag_GetImpScTime = false;
        }
        if (band.use_analytic_band) 
             ImpScGamma = band.GetAnalyticImpurityRate(par_iter->energy, DA, Rho, eps[SILICON], frickel);
        else 
             ImpScGamma = GetImpScGamma();
        
        if (ImpScGamma > 1.0e-20) ImpScTf = imprnl / ImpScGamma; else ImpScTf = 1.0e99;
    }
    else {
        ImpScTf = 2 * dt;
        ImpScGamma = 1.0 / scrt0;
    }

    if((Flag_SurfaceScatter) && ((*c_InSurfRegion)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)])) {
        if(Flag_GetSurfScTime) {
            ssnl = -log(Random());
            Flag_GetSurfScTime = false;
        }
        GetSurfScRate();
        if(SurfScGamma > 0 && ssnl > 0) SurfScTf = ssnl / SurfScGamma; else SurfScTf = 2 * dt;
    } else {
        SurfScTf = 2 * dt;
    }
    
    old_flag = flag;
    Tf = Min(TetTf, CellTf, PhScTf, ImpScTf, SurfScTf, left_time, flag);

    // 记录计算出的时间，用于下一次崩溃时的回溯
    curr_step.calc_TetTf = TetTf;
    curr_step.calc_CellTf = CellTf;
    curr_step.trigger_flag = flag;

    // [DEBUG] 崩溃现场详细尸检
    if(Tf < MY_ZERO && flag != 2){
        cout << "\n================== CRITICAL DIAGNOSIS ==================" << endl;
        cout << "Error: Tf < 0 (Deadlock) Detected at Loop " << loop << endl;
        cout << "old flag" << old_flag << ", new flag " << flag << endl;
        cout << "本轮的速度: (" << vx * velo0 << ", " << vy * velo0 << ", " << vz * velo0 << ") m/s" << endl;
        cout << "本轮的位置: (" << x * spr0 << ", " << y * spr0 << ", " << z * spr0 << ") m" << endl;
        cout << "本轮的 K-矢量: (" << kx * spk0 << ", " << ky * spk0 << ", " << kz * spk0 << ") " << endl;
        cout << "上轮的速度: (" << prev_step.vx * velo0 << ", " << prev_step.vy * velo0 << ", " << prev_step.vz * velo0 << ") m/s" << endl;
        cout << "上轮的位置: (" << prev_step.x * spr0 << ", " << prev_step.y * spr0 << ", " << prev_step.z * spr0 << ") m" << endl;
        cout << "上轮的 K-矢量: (" << prev_step.kx * spk0 << ", " << prev_step.ky * spk0 << ", " << prev_step.kz * spk0 << ") " << endl;
        Flag_Catch = true;
        break;
        
        //cout << "\n[1. REAL SPACE CHECK]" << endl;
        //cout << "Pos Y       : " << y * spr0 << " m" << endl;
        //cout << "Cell Bound L: " << ly[jcell] * spr0 << " m (Index " << jcell << ")" << endl;
        //cout << "Cell Bound R: " << ly[jcell+1] * spr0 << " m" << endl;
        //double dL = y - ly[jcell];
        //double dR = ly[jcell+1] - y;
        //cout << "Dist Left   : " << dL * spr0 << " m" << (fabs(dL) < 1e-12 ? " [ON BOUNDARY]" : "") << endl;
        //cout << "Dist Right  : " << dR * spr0 << " m" << (fabs(dR) < 1e-12 ? " [ON BOUNDARY]" : "") << endl;
        //cout << "Velocity Y  : " << vy * velo0 << " m/s" << endl;
        //
        //cout << "\n[2. K-SPACE CHECK]" << endl;
        //cout << "K Vector    : (" << kx*spk0 << ", " << ky*spk0 << ", " << kz*spk0 << ")" << endl;
        //cout << "TetTf (Time to K-Boundary) : " << TetTf << endl;
        //cout << "CellTf (Time to R-Boundary): " << CellTf << endl;
        //
        //cout << "\n[3. EVOLUTION HISTORY]" << endl;
        //cout << "--- LOOP " << prev_step.loop_idx << " (Previous) ---" << endl;
        //cout << "  Triggered Flag: " << prev_step.trigger_flag << " (1=HitTet, 2=HitCell)" << endl;
        //cout << "  Vy            : " << prev_step.vy * velo0 << endl;
        //cout << "  Ky            : " << prev_step.ky * spk0 << endl;
        //cout << "  Y Position    : " << prev_step.y * spr0 << endl;
        //
        //cout << "--- LOOP " << curr_step.loop_idx << " (Current/Crashed) ---" << endl;
        //cout << "  Target Flag   : " << flag << endl;
        //cout << "  Vy            : " << vy * velo0 << " (Note Sign Change?)" << endl;
        //cout << "  Ky            : " << ky * spk0 << endl;
        //
        //double dKy = curr_step.ky - prev_step.ky;
        //cout << "\n[4. LOGIC TEST]" << endl;
        //cout << "Delta Ky      : " << dKy * spk0 << endl;
        //cout << "Was TetTf small enough? " << (TetTf < CellTf ? "Yes, K-bound closer" : "No, Real-bound closer") << endl;
        //if (flag == 1) cout << "System WANTED to trigger HitTet, but Tf was 0!" << endl;
        //if (flag == 2) cout << "System WANTED to trigger HitCell (Spatial), meaning K-space didn't hit boundary yet." << endl;
//
        //cout << "========================================================" << endl;
        //Flag_Catch = true;
        //break;
    }
      
    // --- 2. 漂移更新 ---
    TetTf -= Tf;
    CellTf -= Tf;
    PhScTf -= Tf;
    ImpScTf -= Tf;
    SurfScTf -= Tf;
    left_time -= Tf;
    
    phrnl -= band.gamtet[itet] * Tf;
    imprnl -= ImpScGamma * Tf;
    ssnl   -= SurfScGamma * Tf;
    
    kx += Ex * Tf;
    ky += Ey * Tf;
    kz += Ez * Tf;
    
    energy += (vx * Ex + vy * Ey + vz * Ez) * Tf;

    x += vx * Tf;
    y += vy * Tf;
    z += vz * Tf;
   
    if(energy < 0) {
      Flag_Catch = true;
      break;
    }
    
    // 更新历史记录
    prev_step = curr_step;
    old_cellflag = band.flag_cellhit;

    // --- 3. 事件处理 ---
    switch (flag) {
    case 1: { // HitTet (K-space boundary)
      if (band.use_analytic_band) {
          HitAnalyticKGrid();
          if (!Flag_Catch) {
              // 同步到当前粒子状态
              kx = par_iter->kx;
              ky = par_iter->ky;
              kz = par_iter->kz;
              energy = par_iter->energy;
              vx = band.analytic_vx;
              vy = band.analytic_vy;
              vz = band.analytic_vz;
              Flag_GetTetTime = true;
              Flag_GetCellTime = true;
          }
      } else {
          HitTet();
          Flag_GetTetTime = true;
          Flag_GetCellTime = true;
      }
      break;
    }
    case 2: { // HitCell (Real-space boundary)
      if (HitCell()) {
        Flag_GetTetTime = true;
        Flag_GetCellTime = true;
      } else {
        trace_back();
        fly_too_far = true;
      }
      break;
    }
    case 3: { // Phonon Scatter
      if (par_type == PELEC) { 
          if (band.use_analytic_band) {
              band.AnalyticPhononScatter(&(*par_iter));
              Flag_SelfScatter = band.analytic_self_scatter;
              // 散射后立即同步
              kx = par_iter->kx; ky = par_iter->ky; kz = par_iter->kz;
              energy = par_iter->energy;
              band.GetAnalyticV_FromTable(&(*par_iter));
              vx = band.analytic_vx; vy = band.analytic_vy; vz = band.analytic_vz;
          } else {
              ElectronPhononScatter();
          }
      } else if (par_type == PHOLE) { 
          HolePhononScatter();
      }

      if (!Flag_SelfScatter) {
          Flag_GetTetTime = true;
          Flag_GetCellTime = true;
          sttt.phononScatter ++;
      }
      Flag_GetPhScTime = true;
      break;
    }
    case 4 : { // Impurity Scatter
      OutPar(&old_par);
      if(par_type == PELEC) {
          if (band.use_analytic_band) {
              band.AnalyticImpurityScatter(&(*par_iter), DA, Rho, eps[SILICON], frickel, ImpScGamma);
              Flag_SelfScatter = band.analytic_self_scatter;
              kx = par_iter->kx; ky = par_iter->ky; kz = par_iter->kz;
              energy = par_iter->energy;
              band.GetAnalyticV_FromTable(&(*par_iter));
              vx = band.analytic_vx; vy = band.analytic_vy; vz = band.analytic_vz;
          } else {
              ElectronImpurityScatter();
          }
      }
      OutPar(&new_par);
      
      old_flag_getTetTime = Flag_SelfScatter;
      if(!Flag_SelfScatter) {
          Flag_GetTetTime = true;
          Flag_GetCellTime = true;
          sttt.impurityScatter ++;
      }
      Flag_GetImpScTime = true;
      break;
    }
    case 5 : { // Surface Scatter
      ParticleSurfaceScatter();
      if(!Flag_SelfScatter) {
          Flag_GetTetTime = true;
          Flag_GetCellTime = true;
      }
      Flag_GetSurfScTime = true;
      break;
    }
   } // end switch

   if(fly_too_far || Flag_Catch) break;
  } // end while

  if(Flag_Catch) {
      par_iter->i = - 9999;
      catch_par_num ++;
  } else {
    OutPar(par_iter);
    if (!fly_too_far) {
      par_iter->left_time = dt;
      par_iter->flag = 1;
    } else 
      par_iter->flag = 0;
  }

  if (loop > 1000) {
    err_message(TOO_MANY_LOOPS, "particle fly");
    dump_par_info();
  }
}

// 解析能带：处理 K 空间网格碰撞
void MeshQuantities::HitAnalyticKGrid() {
    int dir = band.last_k_col_dir;

    switch (dir) {
        case 0: par_iter->kx_idx++; break; // +kx
        case 1: par_iter->kx_idx--; break; // -kx
        case 2: par_iter->ky_idx++; break; // +ky
        case 3: par_iter->ky_idx--; break; // -ky
        case 4: par_iter->kz_idx++; break; // +kz
        case 5: par_iter->kz_idx--; break; // -kz
        default:
            if (mpi_rank == 0) {
                cout << "Error: HitAnalyticKGrid invalid dir " << dir << endl;
            }
            return;
    }

    int N = band.num_ticks_axis;
    bool out_of_bounds = (par_iter->kx_idx < 0 || par_iter->kx_idx >= N ||
                          par_iter->ky_idx < 0 || par_iter->ky_idx >= N ||
                          par_iter->kz_idx < 0 || par_iter->kz_idx >= N);
    if (out_of_bounds) {
        Flag_Catch = true;
        return;
    }

    // 将 k 定位到新的格点
    par_iter->kx = band.k_ticks_code[par_iter->kx_idx];
    par_iter->ky = band.k_ticks_code[par_iter->ky_idx];
    par_iter->kz = band.k_ticks_code[par_iter->kz_idx];

    // 查表更新速度和能量
    band.GetAnalyticV_By_Index(&(*par_iter));
}

// 审计粒子：统计有效区与 Ghost/无效区的数量
void MeshQuantities::Audit_Particles() {
    long long total_in_memory = 0;
    long long ghost_particles = 0;
    long long valid_particles = 0;

    for (size_t idx = 0; idx < par_list.size(); idx++) {
        int count = par_list[idx].size();
        total_in_memory += count;

        if (count > 0) {
            auto it = par_list[idx].begin();
            int p_i = it->i;
            int p_j = it->j;
            int p_k = it->k;

            bool is_valid = (p_i >= c_ibegin && p_i <= c_iend) &&
                            (p_j >= c_jbegin && p_j <= c_jend) &&
                            (p_k >= c_kbegin && p_k <= c_kend);

            if (is_valid) {
                valid_particles += count;
            } else {
                ghost_particles += count;
                if (ghost_particles < 50 || ghost_particles % 10000 == 0) {
                    cout << "[ZOMBIE FOUND] Cell Index: " << idx
                         << " (Par Coord: " << p_i << "," << p_j << "," << p_k << ")"
                         << " Count: " << count << endl;
                }
            }
        }
    }

    cout << "================ PARTICLE AUDIT ================" << endl;
    cout << "Total In Memory : " << total_in_memory << endl;
    cout << "Valid Region    : " << valid_particles << " (Should match 'actual')" << endl;
    cout << "Ghost/Invalid   : " << ghost_particles << " (The MISSING ones!)" << endl;
    cout << "================================================" << endl;
}

/* -------------------------------------------------------------------------- */
/** @brief 读出要模拟粒子的信息.
 * * 
 * @param par_iter 当前要模拟的粒子.
 */
/* ---------------------------------------------------------------------------- */
void MeshQuantities::InPar(list<Particle>::iterator par_iter)
{
  icell   = par_iter->i;
  jcell   = par_iter->j;
  kcell   = par_iter->k;
  itet    = par_iter->itet;
  isym    = par_iter->isym;
  par_type = par_iter->par_type;
  seed     = par_iter->seed;
  par_id  = par_iter->par_id;

  x       = par_iter->x;
  y       = par_iter->y;
  z       = par_iter->z;
  kx      = par_iter->kx;
  ky      = par_iter->ky;
  kz      = par_iter->kz;
  left_time = par_iter->left_time;
  charge  = par_iter->charge;
  energy  = par_iter->energy;
}
/* -------------------------------------------------------------------------- */
/** @brief 将粒子信息保存
 * 
 * @param par_iter 要保存的粒子迭代器 
 */
/* ---------------------------------------------------------------------------- */
  void MeshQuantities::OutPar(Particle * par_iter)
{
  par_iter->i = icell;
  par_iter->j = jcell;
  par_iter->k = kcell;
  par_iter->itet = itet;
  par_iter->isym = isym;
  par_iter->par_type = par_type;
  par_iter->seed = seed;

  par_iter->x = x;
  par_iter->y = y;
  par_iter->z = z;
  par_iter->kx = kx;
  par_iter->ky = ky;
  par_iter->kz = kz;
  par_iter->left_time = left_time;
  par_iter->charge = charge;
  par_iter->energy = energy;
}


  void MeshQuantities::OutPar(list<Particle>::iterator par_iter)
{
  par_iter->i = icell;
  par_iter->j = jcell;
  par_iter->k = kcell;
  par_iter->itet = itet;
  par_iter->isym = isym;
  par_iter->par_type = par_type;
  par_iter->seed = seed;

  par_iter->x = x;
  par_iter->y = y;
  par_iter->z = z;
  par_iter->kx = kx;
  par_iter->ky = ky;
  par_iter->kz = kz;
  par_iter->left_time = left_time;
  par_iter->charge = charge;
  par_iter->energy = energy;
}

void MeshQuantities::compute_field() {

//  if (Flag_compute_potential) {
//if (Flag_NonLinearPoisson)
    linear_poisson_solver();

    save_poisson_pot();

    pot2field_for_hole(p_poisson_pot_saved);

    if ((Flag_QuantumCorrection) && (step >= quantum_start_step)) {
      quantumPotential();
    }
//  }
/*
#ifdef DEBUG_DATA_PRINT
    if (step % debug_print_step == 0) 
    string filename = getFileName("poipot", step);
    print_p_data(p_poisson_pot, filename, pot0);
    }
#endif 
*/
  //pot2field(Flag_QuantumCorrection);
    pot2field(0, p_poisson_pot);
}

void MeshQuantities::print_qc_pot() {

  string filename = getFileName("qpot", step);
  print_p_data(stat_qc_pot, filename, pot0 / quantum_print_steps);
  stat_qc_pot->PutScalar(0);
}

void MeshQuantities::quantum_stat_pot() {
  int i,j,k;
  double *pot, *qstat_pot;

  p_poisson_pot->ExtractView(&pot);
  p_quantum_stat_pot->ExtractView(&qstat_pot);

  for(i = p_ibegin; i <= p_iend; i ++)
    for(k = p_kbegin; k <= p_kend; k ++)
      for(j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap;j ++)
	qstat_pot[P_LINDEX_ONE(i,j,k)] += pot[P_LINDEX_ONE(i,j,k)];
}

void MeshQuantities::clear_quantum_stat_pot() {
  int i,j,k;
  double *qstat_pot;

  p_quantum_stat_pot->ExtractView(&qstat_pot);

  for(i = p_ibegin; i <= p_iend; i ++)
    for(k = p_kbegin; k <= p_kend; k ++)
      for(j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap;j ++)
	qstat_pot[P_LINDEX_ONE(i,j,k)] = 0;
}


void MeshQuantities::save_poisson_pot() {
  int i,j,k;
  double *pot, *saved_pot;

  p_poisson_pot->ExtractView(&pot);
  p_poisson_pot_saved->ExtractView(&saved_pot);

  for(i = p_ibegin; i <= p_iend; i ++)
    for(k = p_kbegin; k <= p_kend; k ++)
      for(j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap;j ++)
	saved_pot[P_LINDEX_ONE(i,j,k)] = pot[P_LINDEX_ONE(i,j,k)];
}

void MeshQuantities::print_subband_energy() {
   int i,j, val;
   MPI_Status status;

   if (mpi_rank != 0)
     MPI_Recv(&j, 1, MPI_INT, mpi_rank - 1, 99 , MPI_COMM_WORLD, &status);

   ofstream ofile;

   string filename = getFileName("subbands", step);

   filename = "./data/" + filename;
   ofile.open(filename.c_str(), iostream::app);

   for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
     for (val = 0; val < valley_num; val ++) 
	 for (i = 0; i < max_subband; i ++) 
	   ofile << j << ' ' << val << ' ' << i << ' ' << subbands[P_SHIFT_NONOVERLAP_Y(j)][val][i] << endl;
     
   ofile.close();

   if (mpi_rank != mpi_size - 1)
     MPI_Send(&j, 1, MPI_INT, mpi_rank + 1, 99, MPI_COMM_WORLD);
}


double MeshQuantities::effective_potential(int pi, int pj, int pk) {
  int i,j,k;
  double ret = 0;
  double x1, x2, y1, y2, z1, z2;
  int * pmat;
  double * poi_pot, tmp_pot;
  double * range_x1, * range_x2, * range_y1, * range_y2, * range_z1, *range_z2;
  double * qc_Eb;

  p_qc_mat->ExtractView(&pmat);
  p_qc_pot->ExtractView(&poi_pot);

  p_qc_range_x1->ExtractView(&range_x1);
  p_qc_range_x2->ExtractView(&range_x2);
  p_qc_range_y1->ExtractView(&range_y1);
  p_qc_range_y2->ExtractView(&range_y2);
  p_qc_range_z1->ExtractView(&range_z1);
  p_qc_range_z2->ExtractView(&range_z2);
  p_qc_Eb->ExtractView(&qc_Eb);

  for (i = lower_x[pi]; i <= upper_x[pi]; i ++)
    for (j = lower_y[pj]; j <= upper_y[pj]; j ++)
      for (k = lower_z[pk]; k <= upper_z[pk]; k ++)
        if ((pmat[P_QCLINDEX(i,j,k)] & NODE_OXIDE) || (pmat[P_QCLINDEX(i,j,k)] & NODE_SILICON)) {

	  tmp_pot = poi_pot[P_QCLINDEX(i,j,k)] + qc_Eb[P_QCLINDEX(i,j,k)];

	  x1 = range_x1[P_QCLINDEX(i,j,k)] - qclx[pi];
	  x2 = range_x2[P_QCLINDEX(i,j,k)] - qclx[pi];
	  y1 = range_y1[P_QCLINDEX(i,j,k)] - qcly[pj];
	  y2 = range_y2[P_QCLINDEX(i,j,k)] - qcly[pj];
	  z1 = range_z1[P_QCLINDEX(i,j,k)] - qclz[pk];
	  z2 = range_z2[P_QCLINDEX(i,j,k)] - qclz[pk];

	  ret += tmp_pot * (erf(z2) - erf(z1)) * (erf(y2) - erf(y1)) * (erf(x2) - erf(x1));
      }
    ret *= 0.125;

  return ret;
}

void MeshQuantities::quantumPotential() {
  int i,j,k;
  double *pot;
  string filename;
  int * pmat;

  p_qc_pot->Import(*p_poisson_pot, *p_qc_Importer, Insert);
  p_poisson_pot->ExtractView(&pot);
  p_qc_mat->ExtractView(&pmat);

  for (j = p_jbegin_nonoverlap;j <= p_jend_nonoverlap; j ++)
//    if (NOT_GHOST_BC(j)){
    {
      for (i = p_ibegin; i <= p_iend; i ++)
        for (k = p_kbegin; k <= p_kend; k ++)
	  if (pmat[P_QCLINDEX(i,j,k)] & NODE_QUANTUM) {
	    pot[P_LINDEX_ONE(i,j,k)] = effective_potential(i,j,k); 
          } 
    }

#ifdef DEBUG_DATA_PRINT
  if (step % debug_print_step == 0) {
    filename = getFileName("ep", step);
    print_p_data(p_poisson_pot, filename, pot0);
    }
#endif
}

void MeshQuantities::init_ep_range() {
  int i,j,k;
  int * pmat;
  double * x1, * x2, * y1, * y2, * z1, * z2;
  double x11, x22, y11, y22, z11, z22;
  double * ep_Eb;

  Epetra_Vector p_range_x1(*p_map_nonoverlap);
  Epetra_Vector p_range_x2(*p_map_nonoverlap);
  Epetra_Vector p_range_y1(*p_map_nonoverlap);
  Epetra_Vector p_range_y2(*p_map_nonoverlap);
  Epetra_Vector p_range_z1(*p_map_nonoverlap);
  Epetra_Vector p_range_z2(*p_map_nonoverlap);
  Epetra_Vector p_Eb(*p_map_nonoverlap);

  p_material->ExtractView(&pmat);
  p_range_x1.ExtractView(&x1);
  p_range_x2.ExtractView(&x2);
  p_range_y1.ExtractView(&y1);
  p_range_y2.ExtractView(&y2);
  p_range_z1.ExtractView(&z1);
  p_range_z2.ExtractView(&z2);
  p_Eb.PutScalar(0);
  p_Eb.ExtractView(&ep_Eb);

  for (j = p_jbegin_nonoverlap;j <= p_jend_nonoverlap; j ++)
//    if (NOT_GHOST_BC(j)){
  {
      for (i = p_ibegin; i <= p_iend; i ++)
        for (k = p_kbegin; k <= p_kend; k ++){
	  if ((pmat[P_LINDEX_ONE(i,j,k)] & NODE_OXIDE) && (!(pmat[P_LINDEX_ONE(i,j,k)] & NODE_SILICON))) 
//              (pmat[P_LINDEX_ONE(i,j,k)] & NODE_UP_BC) || (pmat[P_QCLINDEX(i,j,k)] & NODE_DOWN_BC) ||
//              (pmat[P_LINDEX_ONE(i,j,k)] & NODE_BACK_BC) || (pmat[P_QCLINDEX(i,j,k)] & NODE_FRONT_BC)) 
            ep_Eb[P_LINDEX_ONE(i,j,k)] = -Eb;

	  if (pmat[P_LINDEX_ONE(i,j,k)] & NODE_UP_BC) 
	    x11 = lx[i] - qc_xratio * qc_xtheta;
	  else 
	    x11 = lx[i] - 0.5 * dx[i - 1];
	  if (pmat[P_LINDEX_ONE(i,j,k)] & NODE_DOWN_BC) 
	    x22 = lx[i] + qc_xratio * qc_xtheta;
	  else 
	    x22 = lx[i] + 0.5 * dx[i];
	  if (pmat[P_LINDEX_ONE(i,j,k)] & NODE_LEFT_BC) 
	    y11 = ly[j] - qc_yratio * qc_ytheta;
	  else 
	    y11 = ly[j] - 0.5 * dy[j - 1];
	  if (pmat[P_LINDEX_ONE(i,j,k)] & NODE_RIGHT_BC) 
	    y22 = ly[j] + qc_yratio * qc_ytheta;
	  else 
	    y22 = ly[j] + 0.5 * dy[j];
	  if (pmat[P_LINDEX_ONE(i,j,k)] & NODE_FRONT_BC) 
	    z11 = lz[k] - qc_zratio * qc_ztheta;
	  else 
	    z11 = lz[k] - 0.5 * dz[k - 1];
	  if (pmat[P_LINDEX_ONE(i,j,k)] & NODE_BACK_BC) 
	    z22 = lz[k] + qc_zratio * qc_ztheta;
	  else 
	    z22 = lz[k] + 0.5 * dz[k];

	  x1[P_LINDEX_ONE(i,j,k)] = x11  / (SQRT2 * qc_xtheta);
	  x2[P_LINDEX_ONE(i,j,k)] = x22  / (SQRT2 * qc_xtheta);
	  y1[P_LINDEX_ONE(i,j,k)] = y11  / (SQRT2 * qc_ytheta);
	  y2[P_LINDEX_ONE(i,j,k)] = y22  / (SQRT2 * qc_ytheta);
	  z1[P_LINDEX_ONE(i,j,k)] = z11  / (SQRT2 * qc_ztheta);
	  z2[P_LINDEX_ONE(i,j,k)] = z22  / (SQRT2 * qc_ztheta);
    }
  }

  p_qc_range_x1->Import(p_range_x1, *p_qc_Importer, Insert);
  p_qc_range_x2->Import(p_range_x2, *p_qc_Importer, Insert);
  p_qc_range_y1->Import(p_range_y1, *p_qc_Importer, Insert);
  p_qc_range_y2->Import(p_range_y2, *p_qc_Importer, Insert);
  p_qc_range_z1->Import(p_range_z1, *p_qc_Importer, Insert);
  p_qc_range_z2->Import(p_range_z2, *p_qc_Importer, Insert);
  p_qc_Eb->Import(p_Eb, *p_qc_Importer, Insert);

  qclx.resize(lx.size());
  qcly.resize(ly.size());
  qclz.resize(lz.size());

  for (i = 0;i < lx.size(); i ++)
    qclx[i] = lx[i] / (SQRT2 * qc_xtheta);
  for (i = 0;i < ly.size(); i ++)
    qcly[i] = ly[i] / (SQRT2 * qc_ytheta);
  for (i = 0;i < lz.size(); i ++)
    qclz[i] = lz[i] / (SQRT2 * qc_ztheta);
}

void MeshQuantities::init_p_mat() {
  int i,j,k;
  int * material;
  int * pmat, *qc_mat;
  int * quantumRegion;

  c_material->ExtractView(&material);
  c_quantumRegion->ExtractView(&quantumRegion);

  //p_material->PutScalar(0);
  p_material->ExtractView(&pmat);

  for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
    for (i = 0;i < p_numx; i ++)
      for (k = 0; k < p_numz; k ++) 
        pmat[P_LINDEX_ONE(i,j,k)] = 0;

  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
      for (j = c_jbegin_ghost; j <= c_jend; j ++){

        if (quantumRegion[C_LINDEX_GHOST_ONE(i,j,k)]){
          if (j >= p_jbegin_nonoverlap){
            pmat[P_LINDEX_ONE(i,j,k)] |= NODE_QUANTUM;
            pmat[P_LINDEX_ONE(i + 1,j,k)] |= NODE_QUANTUM;
            pmat[P_LINDEX_ONE(i,j,k + 1)] |= NODE_QUANTUM;
            pmat[P_LINDEX_ONE(i + 1,j,k + 1)] |= NODE_QUANTUM;
          }
          if (j + 1 <= p_jend_nonoverlap){
            pmat[P_LINDEX_ONE(i,j + 1,k)] |= NODE_QUANTUM;
            pmat[P_LINDEX_ONE(i + 1,j + 1,k)] |= NODE_QUANTUM;
            pmat[P_LINDEX_ONE(i,j + 1,k + 1)] |= NODE_QUANTUM;
            pmat[P_LINDEX_ONE(i + 1,j + 1,k + 1)] |= NODE_QUANTUM;
          }
	      }
        if (material[C_LINDEX_GHOST_ONE(i,j,k)] == SILICON){
          if (j >= p_jbegin_nonoverlap){
            pmat[P_LINDEX_ONE(i,j,k)] |= NODE_SILICON;
            pmat[P_LINDEX_ONE(i + 1,j,k)] |= NODE_SILICON;
            pmat[P_LINDEX_ONE(i,j,k + 1)] |= NODE_SILICON;
            pmat[P_LINDEX_ONE(i + 1,j,k + 1)] |= NODE_SILICON;
          }
          if (j + 1 <= p_jend_nonoverlap){
            pmat[P_LINDEX_ONE(i,j + 1,k)] |= NODE_SILICON;
            pmat[P_LINDEX_ONE(i + 1,j + 1,k)] |= NODE_SILICON;
            pmat[P_LINDEX_ONE(i,j + 1,k + 1)] |= NODE_SILICON;
            pmat[P_LINDEX_ONE(i + 1,j + 1,k + 1)] |= NODE_SILICON;
          }
	      }
        if (material[C_LINDEX_GHOST_ONE(i,j,k)] == OXIDE){
          if (j >= p_jbegin_nonoverlap){
            pmat[P_LINDEX_ONE(i,j,k)] |= NODE_OXIDE;
            pmat[P_LINDEX_ONE(i + 1,j,k)] |= NODE_OXIDE;
            pmat[P_LINDEX_ONE(i,j,k + 1)] |= NODE_OXIDE;
            pmat[P_LINDEX_ONE(i + 1,j,k + 1)] |= NODE_OXIDE;
          }
          if (j + 1 <= p_jend_nonoverlap){
            pmat[P_LINDEX_ONE(i,j + 1,k)] |= NODE_OXIDE;
            pmat[P_LINDEX_ONE(i + 1,j + 1,k)] |= NODE_OXIDE;
            pmat[P_LINDEX_ONE(i,j + 1,k + 1)] |= NODE_OXIDE;
            pmat[P_LINDEX_ONE(i + 1,j + 1,k + 1)] |= NODE_OXIDE;
          }
	      }
      }

  p_qc_mat->Import(*p_material, *p_qc_Importer, Insert);
  p_qc_mat->ExtractView(&qc_mat);

  for (i = 0;i < p_numx; i ++)
    for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      for (k = 0; k < p_numz; k ++) {
        if ((i == 0) || (!((qc_mat[P_QCLINDEX(i - 1, j , k)] & NODE_OXIDE) || (qc_mat[P_QCLINDEX(i - 1, j , k)] & NODE_SILICON))))
          pmat[P_LINDEX_ONE(i,j,k)] |= NODE_UP_BC;
        if ((i == p_numx - 1) || (!((qc_mat[P_QCLINDEX(i + 1, j , k)] & NODE_OXIDE) || (qc_mat[P_QCLINDEX(i + 1, j , k)] & NODE_SILICON))))
          pmat[P_LINDEX_ONE(i,j,k)] |= NODE_DOWN_BC;
        if ((k == 0) || (!((qc_mat[P_QCLINDEX(i , j , k - 1)] & NODE_OXIDE) || (qc_mat[P_QCLINDEX(i , j , k - 1)] & NODE_SILICON))))
          pmat[P_LINDEX_ONE(i,j,k)] |= NODE_FRONT_BC;
        if ((k == p_numz - 1) || (!((qc_mat[P_QCLINDEX(i , j , k + 1)] & NODE_OXIDE) || (qc_mat[P_QCLINDEX(i , j , k + 1)] & NODE_SILICON))))
          pmat[P_LINDEX_ONE(i,j,k)] |= NODE_BACK_BC;
        if ((j == 0) || (!((qc_mat[P_QCLINDEX(i , j - 1, k)] & NODE_OXIDE) || (qc_mat[P_QCLINDEX(i , j - 1, k )] & NODE_SILICON))))
                  pmat[P_LINDEX_ONE(i,j,k)] |= NODE_LEFT_BC;
        if ((j == p_numy - 1) || (!((qc_mat[P_QCLINDEX(i , j + 1 , k)] & NODE_OXIDE) || (qc_mat[P_QCLINDEX(i , j + 1, k)] & NODE_SILICON))))
          pmat[P_LINDEX_ONE(i,j,k)] |= NODE_RIGHT_BC;
      }

  p_qc_mat->Import(*p_material, *p_qc_Importer, Insert);

}

void MeshQuantities::init_effective_potential() {
  int i,j,k, sub;
  int * p_global_id;
  vector<double>::iterator low_iter;
  vector<double>::iterator uper_iter;

  low_iter = lx.begin();
  uper_iter = lx.begin();

  for (i = 0; i < p_numx; i ++){
    low_iter = lower_bound(low_iter, lx.begin() + i, lx[i] - qc_xratio * qc_xtheta);
    uper_iter = upper_bound(uper_iter, lx.end(), lx[i] + qc_xratio * qc_xtheta);

    if (low_iter != lx.begin()) low_iter --;
    if (uper_iter == lx.end()) uper_iter --;

    lower_x.push_back(low_iter - lx.begin());
    upper_x.push_back(uper_iter - lx.begin());
  }

  low_iter = ly.begin();
  uper_iter = ly.begin();

  for (i = 0; i < p_numy; i ++){
    low_iter = lower_bound(low_iter, ly.begin() + i, ly[i] - qc_yratio * qc_ytheta);
    uper_iter = upper_bound(uper_iter, ly.end(), ly[i] + qc_yratio * qc_ytheta);

    if (low_iter != ly.begin()) low_iter --;
    if (uper_iter == ly.end()) uper_iter --;

    lower_y.push_back(low_iter - ly.begin());
    upper_y.push_back(uper_iter - ly.begin());
  }

  low_iter = lz.begin();
  uper_iter = lz.begin();

  for (i = 0; i < p_numz; i ++){
    low_iter = lower_bound(low_iter, lz.begin() + i, lz[i] - qc_zratio * qc_ztheta);
    uper_iter = upper_bound(uper_iter, lz.end(), lz[i] + qc_zratio * qc_ztheta);

    if (low_iter != lz.begin()) low_iter --;
    if (uper_iter == lz.end()) uper_iter --;

    lower_z.push_back(low_iter - lz.begin());
    upper_z.push_back(uper_iter - lz.begin());
  }

  p_qc_jbegin = lower_y[p_jbegin_nonoverlap];
  p_qc_jend = upper_y[p_jend_nonoverlap];
  p_qc_numy = (p_qc_jend - p_qc_jbegin + 1) * p_numxz;
  

//  cout << "mpi_rank = " << p_rank << ' ' << p_qc_jbegin << ' ' << p_qc_jend << ' ' << p_qc_numy << endl;

  //cout << "mpi_rank = " << mpi_rank << ' ' << p_qc_jbegin <<  ' ' << p_jbegin_nonoverlap << ' ' << p_qc_jend << ' ' << p_jend_nonoverlap << endl;

  p_global_id = (int *) malloc(sizeof(int) * p_qc_numy);
  
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_qc_jbegin; j <= p_qc_jend; j ++)
	p_global_id[P_QCLINDEX(i,j,k)] = P_GINDEX(i,j,k);

  p_qc_map = new Epetra_Map(-1, p_qc_numy, p_global_id, 0, *Comm);
  p_qc_mat = new Epetra_IntVector(*p_qc_map);
  p_qc_pot = new Epetra_Vector(*p_qc_map);

  p_qc_range_x1 = new Epetra_Vector(*p_qc_map);
  p_qc_range_x2 = new Epetra_Vector(*p_qc_map);
  p_qc_range_y1 = new Epetra_Vector(*p_qc_map);
  p_qc_range_y2 = new Epetra_Vector(*p_qc_map);
  p_qc_range_z1 = new Epetra_Vector(*p_qc_map);
  p_qc_range_z2 = new Epetra_Vector(*p_qc_map);
  p_qc_Eb = new Epetra_Vector(*p_qc_map);

  p_qc_Importer = new Epetra_Import(*p_qc_map, *p_map_nonoverlap);

}

void MeshQuantities::compute_cell_charge() {

   double * c_electron_charge_value, * c_hole_charge_value, *cvol, *c_da_value, *c_par_charge_value;
   int i,j,k;
   double charge;
   c_volume->ExtractView(&cvol);
   c_da->ExtractView(&c_da_value);

   list<Particle> * c_par_list;

   list<Particle>::iterator iter;

   c_electron_charge->PutScalar(0);
   c_hole_charge->PutScalar(0);
   c_par_charge->PutScalar(0);

   c_electron_charge->ExtractView(&c_electron_charge_value);
   c_hole_charge->ExtractView(&c_hole_charge_value);
   c_par_charge->ExtractView(&c_par_charge_value);

    /* loop for each entry */
   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
      /*including the ghost cells */
      for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++)
      { 
	      /*list of particles in cell (i,j) */
        c_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];

        /*loop for each particles in this cell*/
        for (iter = c_par_list->begin(); iter != c_par_list->end(); iter ++){
          charge = iter->charge;
 
	        if (iter->par_type == PELEC)
	          c_electron_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] += charge;
	        else 
	          c_hole_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] += charge;

	        c_par_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] += fabs(charge);
       }
      }

     for (i = c_ibegin; i <= c_iend; i ++)
        for (k = c_kbegin; k <= c_kend; k ++)
          for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++){
	          c_par_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] /= cvol[C_LINDEX_GHOST_ONE(i,j,k)];
	          c_par_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] = Max(c_par_charge_value[C_LINDEX_GHOST_ONE(i,j,k)], 0.5 * c_da_value[C_LINDEX_GHOST_ONE(i,j,k)]);
          }
}

void MeshQuantities::linear_poisson_solver() {

  p_poisson_pot->PutScalar(0);

  compute_rhs();

  t_time->ResetStartTime();

  poisson_solver->solve_poisson(p_rhs, p_poisson_pot);
   //you can plot and check the poisson results 

  if (mpi_rank == 0) {
    cout << "time used: " << t_time->ElapsedTime() << "s" << endl;
    //  << DSolver.NumIters() 
//         << " steps, residual's norm = " << DSolver.TrueResidual() << endl
  }
 }

 void MeshQuantities::pot2field_for_hole(Epetra_Vector * p_pot_vec) {
   int i,j,k;
   double *field_x, *field_y,*field_z, *pot, *qc_pot;

   /* we may not have enough potential value need to compute field value, so
    * import it from the neibough process */ 
   p_pot->Import(*p_pot_vec, *p_Importer, Insert);

   c_work1->ExtractView(&field_x);
   c_work2->ExtractView(&field_y);
   c_work3->ExtractView(&field_z);

   p_pot->ExtractView(&pot);

   /* compute the field value for each cell on the local proceess */
   // page 5 in basic.pdf
   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin; j <= c_jend; j ++)
       {
        field_x[C_LINDEX_ONE(i,j,k)] = -(pot[P_LINDEX_ONE(i + 1,j,k)] - pot[P_LINDEX_ONE(i, j ,k)]
                + pot[P_LINDEX_ONE(i + 1,j+1,k)] - pot[P_LINDEX_ONE(i, j + 1,k)] 
                + pot[P_LINDEX_ONE(i + 1,j,k+1)] - pot[P_LINDEX_ONE(i, j ,k+1)]
                + pot[P_LINDEX_ONE(i + 1,j+1,k+1)] - pot[P_LINDEX_ONE(i, j + 1,k+1)]) / (4.0 * dx[i]);

        field_y[C_LINDEX_ONE(i,j,k)] = -(pot[P_LINDEX_ONE(i ,j + 1,k)] + pot[P_LINDEX_ONE(i + 1, j + 1,k)]
                +pot[P_LINDEX_ONE(i ,j + 1,k+1)] + pot[P_LINDEX_ONE(i + 1, j + 1,k+1)]
                - pot[P_LINDEX_ONE(i ,j,k)] - pot[P_LINDEX_ONE(i + 1, j,k)]
                - pot[P_LINDEX_ONE(i ,j,k+1)] - pot[P_LINDEX_ONE(i + 1, j,k+1)]) / (4.0 * dy[j]);
        field_z[C_LINDEX_ONE(i,j,k)] = -(pot[P_LINDEX_ONE(i ,j,k+1)] + pot[P_LINDEX_ONE(i + 1, j ,k+1)]
                +pot[P_LINDEX_ONE(i ,j + 1,k+1)] + pot[P_LINDEX_ONE(i + 1, j + 1,k+1)]
                - pot[P_LINDEX_ONE(i ,j,k)] - pot[P_LINDEX_ONE(i + 1, j,k)]
                - pot[P_LINDEX_ONE(i ,j+1,k)] - pot[P_LINDEX_ONE(i + 1, j+1,k)]) / (4.0 * dz[k]);
     }

   /* field values in the ghost cells are aslo needed while simulating
    * particles , import them */
   c_h_field_x->Import(*c_work1, *c_Importer, Insert);
   c_h_field_y->Import(*c_work2, *c_Importer, Insert);
   c_h_field_z->Import(*c_work3, *c_Importer, Insert);

 }


 void MeshQuantities::pot2field(int flag, Epetra_Vector * p_pot_vec) {
   int i,j,k;
   double *field_x, *field_y,*field_z, *pot, *qc_pot;
   int * quantumRegion;

   c_quantumRegion->ExtractView(&quantumRegion);

   /* we may not have enough potential value need to compute field value, so
    * import it from the neibough process */ 
   p_pot->Import(*p_pot_vec, *p_Importer, Insert);

   c_work1->ExtractView(&field_x);
   c_work2->ExtractView(&field_y);
   c_work3->ExtractView(&field_z);

   p_pot->ExtractView(&pot);

   /*compute the field value for each cell on the local proceess*/
   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin; j <= c_jend; j ++)
       if ((!flag) || (quantumRegion[C_LINDEX_ONE(i,j,k)]))
//       if (NOT_GHOST_CELL(j))
       {
	 field_x[C_LINDEX_ONE(i,j,k)] = -(pot[P_LINDEX_ONE(i + 1,j,k)] - pot[P_LINDEX_ONE(i, j ,k)]
					+ pot[P_LINDEX_ONE(i + 1,j+1,k)] - pot[P_LINDEX_ONE(i, j + 1,k)] 
					+ pot[P_LINDEX_ONE(i + 1,j,k+1)] - pot[P_LINDEX_ONE(i, j ,k+1)]
					+ pot[P_LINDEX_ONE(i + 1,j+1,k+1)] - pot[P_LINDEX_ONE(i, j + 1,k+1)]) / (4.0 * dx[i]);

	 field_y[C_LINDEX_ONE(i,j,k)] = -(pot[P_LINDEX_ONE(i ,j + 1,k)] + pot[P_LINDEX_ONE(i + 1, j + 1,k)]
					+pot[P_LINDEX_ONE(i ,j + 1,k+1)] + pot[P_LINDEX_ONE(i + 1, j + 1,k+1)]
					- pot[P_LINDEX_ONE(i ,j,k)] - pot[P_LINDEX_ONE(i + 1, j,k)]
					- pot[P_LINDEX_ONE(i ,j,k+1)] - pot[P_LINDEX_ONE(i + 1, j,k+1)]) / (4.0 * dy[j]);
	 field_z[C_LINDEX_ONE(i,j,k)] = -(pot[P_LINDEX_ONE(i ,j,k+1)] + pot[P_LINDEX_ONE(i + 1, j ,k+1)]
					+pot[P_LINDEX_ONE(i ,j + 1,k+1)] + pot[P_LINDEX_ONE(i + 1, j + 1,k+1)]
					- pot[P_LINDEX_ONE(i ,j,k)] - pot[P_LINDEX_ONE(i + 1, j,k)]
					- pot[P_LINDEX_ONE(i ,j+1,k)] - pot[P_LINDEX_ONE(i + 1, j+1,k)]) / (4.0 * dz[k]);
     }

   /*
   if (c_jbegin == 0) 
     for (i = c_ibegin; i <= c_iend; i ++)
       for (k = c_kbegin; k <= c_kend; k ++){
	 field_x[C_LINDEX_ONE(i,0,k)] = field_x[C_LINDEX_ONE(i,1, k)];
	 field_y[C_LINDEX_ONE(i,0,k)] = field_y[C_LINDEX_ONE(i,1, k)];
	 field_z[C_LINDEX_ONE(i,0,k)] = field_z[C_LINDEX_ONE(i,1, k)];
     }

   if (c_jend == c_numy - 1) 
     for (i = c_ibegin; i <= c_iend; i ++)
       for (k = c_kbegin; k <= c_kend; k ++){
	 field_x[C_LINDEX_ONE(i,c_jend,k)] = field_x[C_LINDEX_ONE(i,c_jend - 1, k)];
	 field_y[C_LINDEX_ONE(i,c_jend,k)] = field_y[C_LINDEX_ONE(i,c_jend - 1, k)];
	 field_z[C_LINDEX_ONE(i,c_jend,k)] = field_z[C_LINDEX_ONE(i,c_jend - 1, k)];
       }
       */
   /* field values in the ghost cells are aslo needed while simulating
    * particles , import them */
   c_field_x->Import(*c_work1, *c_Importer, Insert);
   c_field_y->Import(*c_work2, *c_Importer, Insert);
   c_field_z->Import(*c_work3, *c_Importer, Insert);

 #ifdef DEBUG_DATA_PRINT
/*
  if (step % debug_print_step == 0) {
   string filename;
   filename = getFileName("totpot", step);
   print_p_data(p_poisson_pot, filename, pot0);
   string filename;

   filename = getFileName("xfield", step);
   print_c_data(c_field_x, filename, false, field0);
   filename = getFileName("yfield", step);
   print_c_data(c_field_y, filename, false, field0);
   filename = getFileName("zfield", step);
   print_c_data(c_field_z, filename, false, field0);
  }
*/
 #endif
 }

/*compute the rhs for poisson equation*/
 void MeshQuantities::compute_rhs() {

   double * p_par_charge_value, * p_rhs_value, * p_surf_value, * p_dop_value;
   double *vol_val;
   int i,j,k,icont, iplane;
   double * vadd_val;
   
   // 提取电势与材料属性
   double * pot_val;
   int * p_mat_val;

   p_vadd->ExtractView(&vadd_val);
   p_dop->ExtractView(&p_dop_value);
   p_volume->ExtractView(&vol_val);
   
   /*initialize with zero*/
   p_rhs->PutScalar(0);

   p_par_charge->ExtractView(&p_par_charge_value);
   p_rhs->ExtractView(&p_rhs_value);
   p_poisson_pot->ExtractView(&pot_val);
   p_material->ExtractView(&p_mat_val);

   /*only consider the silicon part */
   for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      {
         int idx = P_LINDEX_ONE(i,j,k);
         
         if (band.use_analytic_band) {
             // 解析能带：显式加入空穴
             double rhs_val = p_dop_value[idx] + p_par_charge_value[idx] * vol_val[idx];
             
             if (p_mat_val[idx] & NODE_SILICON) {
                 double phi = pot_val[idx];
                 double p_dens = band.Ni * exp(-phi);
                 double Q_hole = p_dens * vol_val[idx];
                 rhs_val += Q_hole;
             }
             
             p_rhs_value[idx] = rhs_val;
         } else {
             // 全能带：保持原逻辑
             if (fabs(p_par_charge_value[idx]) > 0){
                 p_rhs_value[idx] = 
                   p_par_charge_value[idx] * vol_val[idx] + p_dop_value[idx];
             }
         }
      }

   for(icont=0;icont<contact.size();icont++)
    for(iplane=0;iplane<contact[icont].NumContactPlane;iplane++)
      for(i=contact[icont].BeginI[iplane];i<=contact[icont].EndI[iplane];i++)
    for(k=contact[icont].BeginK[iplane];k<=contact[icont].EndK[iplane];k++)
      for(j=contact[icont].BeginJ[iplane];j<=contact[icont].EndJ[iplane];j++)
        if ((j >= p_jbegin_nonoverlap) && (j <= p_jend_nonoverlap)){
          p_rhs_value[P_LINDEX_ONE(i,j,k)] = contact[icont].CurrentVapp + vadd_val[P_LINDEX_ONE(i,j,k)] + contact[icont].PhiMS;
        }
 }

 void MeshQuantities::init_nonlinear_poisson() {
 }

 void MeshQuantities::init_poisson_matrix() {

   poisson_solver = new Amesos_PoissonSolver(this); 

/*
   int LevelFill = 0, Overlap = 2;

   Ifpack_IlukGraph Graph(A->Graph(), LevelFill, Overlap);

   Graph.ConstructFilledGraph();

   RILU = new Ifpack_CrsRiluk(Graph);

   int initerr = RILU->InitValues(*A);

   RILU->Factor();
*/

 }

 void MeshQuantities::print_particle_data(string filename) {

   MPI_Status status;
   int i,j,k;
   list<Particle> * c_par_list;
   list<Particle>::iterator iter;

   if (mpi_rank != 0)
     MPI_Recv(&j, 1, MPI_INT, mpi_rank - 1, 99 , MPI_COMM_WORLD, &status);

   ofstream ofile;
   filename = "./data/" + filename;
   ofile.open(filename.c_str(), iostream::app);

   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
     for (j = c_jbegin; j <= c_jend; j ++){

       c_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];

       for (iter = c_par_list->begin(); iter != c_par_list->end(); iter ++){
         ofile << iter->par_id << ' ' << endl;
         /*              << iter->i << ' ' << iter->j << ' '
               << iter->x << ' ' << iter->y << ' '
               << iter->kx << ' ' << iter->ky << ' ' << iter->kz << ' '
               << iter->charge << ' ' << iter->energy << ' ' << iter->par_type << endl;*/
       }
     }

   ofile.close();

   if (mpi_rank != mpi_size - 1)
     MPI_Send(&j, 1, MPI_INT, mpi_rank + 1, 99, MPI_COMM_WORLD);
 }

void MeshQuantities::select_kstate_fermi_dirac(Particle * iter, int dir, double Ef) {
  double ddos,rdos;
  int iptype; 
         /*particle type*/
           iptype = iter->par_type;
           seed = iter->seed;
           FindK:
           /*AcRejection to find particle energy
             choose an energy according to Boltzmann distribution
             (ignore density of states)*/
      //     energy = random_fermi_fun(Ef, band.emax);

           energy = random_fermi_fun(Ef, band.emax);

	     //get DOS for [article energy
	     ddos = band.CALDOSSUM(energy, iptype);

	     /*AcRejection step to consider DOS
	       maximum DOS times random number*/
	     rdos = Random() * band.DOSMAX[iptype];

	     /*choose again */
	     if (rdos > ddos) goto FindK;
   
	     for(iband = band.bandof[iptype]; iband < band.bandof[iptype] + band.nband[iptype]; iband++) 
	       {
		 rdos -= band.CALDOS(energy , iband);
		 if (rdos < 0) break;		
	       }

	   isym = ((int)(Random() * 48)) % 48;

           /*get a k-vector with unform probability on equi energy surface*/
           get_state();

           GetV();

	   if (((dir == 1) && (vy < 0)) || ((dir == -1) && (vy > 0)))
             goto FindK;

           if (Flag_SelfScatter) goto FindK;
           /*if no state is found, try again*/
           /*save state */
           iter->itet = itet;
           iter->isym = isym;
           iter->seed = seed;
           iter->kx = kx;
           iter->ky = ky;
           iter->kz = kz;
           iter->energy = energy;
}


void MeshQuantities::select_kstate(Particle * iter, int dir) {
  // 解析能带分支
  if (band.use_analytic_band) {
    double tmp_const = 1 - exp(- band.emax);
    iter->energy = -log(1 - Random() * tmp_const);

    band.SelectAnalyticKState(iter, iter->energy);
    band.GetAnalyticV_FromTable(iter);

    if (((dir == 1) && (band.analytic_vy < 0)) || ((dir == -1) && (band.analytic_vy > 0))) {
      select_kstate(iter, dir);
      return;
    }
    iter->itet = 0;
    return;
  }

  double ddos,rdos;
  double tmp_const = 1 - exp(- band.emax);
  int iptype; 
         /*particle type*/
           iptype = iter->par_type;
           seed = iter->seed;
           FindK:
           /*AcRejection to find particle energy
             choose an energy according to Boltzmann distribution
             (ignore density of states)*/
           energy = -log( 1 - Random() * tmp_const);

           //get DOS for [article energy
           ddos = band.CALDOSSUM(energy, iptype);

           /*AcRejection step to consider DOS
             maximum DOS times random number*/
           rdos = Random() * band.DOSMAX[iptype];

           /*choose again */
           if (rdos > ddos) goto FindK;
           /*choose band index */
           for(iband = band.bandof[iptype]; iband < band.bandof[iptype] + band.nband[iptype]; iband++) 
             {
               rdos -= band.CALDOS(energy , iband);
               if (rdos < 0) break;		
             }

	   isym = ((int)(Random() * 48)) % 48;

           /*get a k-vector with unform probability on equi energy surface*/
           get_state();


           if (Flag_SelfScatter) goto FindK;

           GetV();

	   if (((dir == 1) && (vy < 0)) || ((dir == -1) && (vy > 0)))
             goto FindK;

           /*if no state is found, try again*/
           /*save state */
           iter->itet = itet;
           iter->isym = isym;
           iter->seed = seed;
           iter->kx = kx;
           iter->ky = ky;
           iter->kz = kz;
           iter->energy = energy;
}


#include "particle_init.cpp"  // 粒子初始化相关的实现单独放入文件，便于阅读

/*init the vectors according to device file, 
 * such as dopping density for donor and acceptor, 
 * particles' number */
 void MeshQuantities::init_by_user_input() {

    int i;

    // cmd_list.size 代表着 ldg.txt 中命令的行数，因此下面的for循环代表遍历全部输入命令
    for (i = 0;i < cmd_list.size(); i ++) {
      switch (cmd_list[i].type) { // type 对应 read_device_file 函数中赋给 map op 的值 
      case 1: // 对应 donor 命令, 用 range[0]~range[5] 对应了空间的范围，dbl_param 代表了输入的掺杂浓度，会将其赋值给 c_donor
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].dbl_param[0], c_donor);
        break;
      case 2: // 对应 acceptor
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].dbl_param[0], c_acceptor);
        break;
      case 3: // Corresponds to "region"
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[0], c_material);

        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[1], c_desired_electron_number);

        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[2], c_desired_hole_number);
        break;
      case 4: // Corresponds to "motioncube"
        SetMRCube(cmd_list[i].range[0], cmd_list[i].range[1],
                  cmd_list[i].range[2], cmd_list[i].range[3],
                  cmd_list[i].range[4], cmd_list[i].range[5],         // specify the range of the cube
                  cmd_list[i].int_param[0], cmd_list[i].int_param[1], 
                  cmd_list[i].int_param[2], cmd_list[i].int_param[3],
                  cmd_list[i].int_param[4], cmd_list[i].int_param[5]  // Specify the corresponding motion at each boundary of the cube
        );
        break;
      case 5: // Corresponds to "attachcontact"
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[0], c_attached_contact);

        break;
      case 7: // Corresponds to "motionplane"
        SetMRPlane(cmd_list[i].range[0], cmd_list[i].range[1],
                    cmd_list[i].range[2], cmd_list[i].range[3],
        cmd_list[i].range[4], cmd_list[i].range[5],
                    cmd_list[i].int_param[0], cmd_list[i].int_param[1]);
        break;
      case 8: // Correspond to "parnumber"
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[0], c_desired_electron_number);

        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[1], c_desired_hole_number);

        break;
      case 12:  // Correspond to "ScatterArea"
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          cmd_list[i].int_param[0], c_scatter_type);
        break;
      case 17:  // Correspond to "surface_scatter_range"
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          1, c_InSurfRegion);
      case 18:  // Correspond to "quantumRegion"
        init_vector_entry(cmd_list[i].range[0], cmd_list[i].range[1],
                          cmd_list[i].range[2], cmd_list[i].range[3],
                          cmd_list[i].range[4], cmd_list[i].range[5],
                          1, c_quantumRegion);

     }
   }
 }

 void MeshQuantities::init_surface_roughness() {
    int i,j,k, isurf;
    double dist;
    int *nearestSurf, *EeffDirection;

    c_nearestSurf->ExtractView(&nearestSurf);
    c_EeffDirection->ExtractView(&EeffDirection);

    // 找到离各个方向的 Si/Oxide界面 最近的面，这个面可以是由 ghost cell 组成
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
        for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++){
            dist = 1e99;
            for (isurf = 0; isurf < NumSurface; isurf ++) {
              if ((SurfaceType[isurf] == 0) && (dist > fabs(lx[i] - SurfacePosition[isurf]))) {
                dist = fabs(lx[i] - SurfacePosition[isurf]);
                nearestSurf[C_LINDEX_GHOST_ONE(i,j,k)] = isurf;
              }
              if ((SurfaceType[isurf] == 1) && (dist > fabs(ly[j] - SurfacePosition[isurf]))) {
                dist = fabs(ly[j] - SurfacePosition[isurf]);
                nearestSurf[C_LINDEX_GHOST_ONE(i,j,k)] = isurf;
              }
              if ((SurfaceType[isurf] == 2) && (dist > fabs(lz[k] - SurfacePosition[isurf]))) {
                dist = fabs(lz[k] - SurfacePosition[isurf]);
                nearestSurf[C_LINDEX_GHOST_ONE(i,j,k)] = isurf;
              }
            }

            EeffDirection[C_LINDEX_GHOST_ONE(i,j,k)] = SurfaceType[nearestSurf[C_LINDEX_GHOST_ONE(i,j,k)]];
        }
}

void MeshQuantities::init_cell_data() {

  int i,j,k;

  double * volume_value = NULL;

  /**
   * @details 
   * ExtractView is used to provide a direct pointer to the internal data 
   *  of an object so that you can read or modify the data efficiently. 
   */
  // Gets a pointer to the volume data
  c_volume->ExtractView(&volume_value);

   double * electron_charge, * hole_charge ;
   int  * material;
   double dope;
   double * donor, * acceptor;
   double *da_value;

  // Gets a pointer to the donor doping concentration data
  c_donor->ExtractView(&donor);

  // Gets a pointer to the acceptor doping concentration data
  c_acceptor->ExtractView(&acceptor);

  // 
  c_material->ExtractView(&material);

  c_init_electron_charge->ExtractView(&electron_charge);

  c_init_hole_charge->ExtractView(&hole_charge);

  c_da->ExtractView(&da_value);

   sum_charge_proc[PELEC] = 0;
   sum_charge_proc[PHOLE] = 0;


   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++){
	      volume_value[C_LINDEX_GHOST_ONE(i,j,k)] = dx[i] * dy[j] * dz[k];
        da_value[C_LINDEX_GHOST_ONE(i,j,k)] = donor[C_LINDEX_GHOST_ONE(i,j,k)] + acceptor[C_LINDEX_GHOST_ONE(i,j,k)];
     }

   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++)
       if (material[C_LINDEX_GHOST_ONE(i,j,k)] == SILICON){

         dope = donor[C_LINDEX_GHOST_ONE(i,j,k)] - acceptor[C_LINDEX_GHOST_ONE(i,j,k)];

         if (dope > 0) {

           electron_charge[C_LINDEX_GHOST_ONE(i,j,k)] = -dope * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];

           hole_charge[C_LINDEX_GHOST_ONE(i,j,k)] =  (band.Ni * band.Ni / dope) * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];

         } else if (dope < 0) {

           electron_charge[C_LINDEX_GHOST_ONE(i,j,k)] =  (band.Ni * band.Ni / dope) * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];

           hole_charge[C_LINDEX_GHOST_ONE(i,j,k)] = -dope * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];

         }else {

           electron_charge[C_LINDEX_GHOST_ONE(i,j,k)] = -band.Ni * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];

           hole_charge[C_LINDEX_GHOST_ONE(i,j,k)] =  band.Ni * volume_value[C_LINDEX_GHOST_ONE(i,j,k)];
         }
         if ((j >= c_jbegin) && (j <= c_jend)){
           sum_charge_proc[PELEC] += electron_charge[C_LINDEX_GHOST_ONE(i,j,k)];
           sum_charge_proc[PHOLE] += hole_charge[C_LINDEX_GHOST_ONE(i,j,k)];
         }
       }

   MPI_Allreduce(sum_charge_proc, sum_charge, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 #ifdef DEBUG
   if (mpi_rank == 0)
     cout << sum_charge[PELEC] << ' ' << sum_charge[PHOLE] << endl;
 #endif 
 }

 void MeshQuantities::init_point_data() {

    int i,j,k;

    double * dop, *vol, * donor, * acceptor, *cvol;
    double vol_tmp, da_tmp;
    int * material;
    double *charge_fac;
    double alfa, ndop;
    double * vadd_val;

    p_vadd->ExtractView(&vadd_val);
    
    c_material->ExtractView(&material);

    c_volume->ExtractView(&cvol);

    p_volume->ExtractView(&vol);

    p_dop->ExtractView(&dop);

    p_charge_fac->ExtractView(&charge_fac);
  
    c_donor->ExtractView(&donor);

    c_acceptor->ExtractView(&acceptor);

    p_vadd->PutScalar(0);

    for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
        for (j = p_jbegin; j <= p_jend; j ++){
          dop[P_LINDEX_ONE(i,j,k)] = 0;
          vol[P_LINDEX_ONE(i,j,k)] = 0;
      }
          
    for (i = c_ibegin; i <= c_iend; i ++)
      for (k = c_kbegin; k <= c_kend; k ++)
        for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++)
        
          if (material[C_LINDEX_GHOST_ONE(i,j,k)] == SILICON) {

            vol_tmp = cvol[C_LINDEX_GHOST_ONE(i,j,k)] * 0.125;
            
            da_tmp =  vol_tmp * (donor[C_LINDEX_GHOST_ONE(i,j,k)] - acceptor[C_LINDEX_GHOST_ONE(i,j,k)]);

            //        if (NOT_GHOST_CELL(j)){
            if (j >= p_jbegin) {
              dop[P_LINDEX_ONE(i + 1, j, k)] += da_tmp;
              dop[P_LINDEX_ONE(i, j, k)] += da_tmp;
              dop[P_LINDEX_ONE(i + 1, j, k + 1)] += da_tmp;
              dop[P_LINDEX_ONE(i, j, k + 1)] += da_tmp;
            }
            if (j < p_jend) {
              dop[P_LINDEX_ONE(i + 1,j + 1, k)]  += da_tmp;
              dop[P_LINDEX_ONE(i,j + 1, k)] += da_tmp;
              dop[P_LINDEX_ONE(i + 1,j + 1, k + 1)]  += da_tmp;
              dop[P_LINDEX_ONE(i,j + 1, k + 1)] += da_tmp;
            }
            if (j >= p_jbegin) {
              vol[P_LINDEX_ONE(i + 1,j, k)] += vol_tmp;
              vol[P_LINDEX_ONE(i, j, k)] += vol_tmp;
              vol[P_LINDEX_ONE(i + 1,j, k + 1)] += vol_tmp;
              vol[P_LINDEX_ONE(i, j, k + 1)] += vol_tmp;
            }
            if (j < p_jend) {
              vol[P_LINDEX_ONE(i + 1,j + 1, k)]  += vol_tmp;
              vol[P_LINDEX_ONE(i,j + 1, k)] += vol_tmp;
              vol[P_LINDEX_ONE(i + 1,j + 1, k+1)]  += vol_tmp;
              vol[P_LINDEX_ONE(i,j + 1, k + 1)] += vol_tmp;
            }
    //}
        }
  
      for (i = p_ibegin; i <= p_iend; i ++)
        for (k = p_kbegin; k <= p_kend; k ++)
          for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++) {

            if(fabs(vol[P_LINDEX_ONE(i,j,k)]) > 1e-16 )
              ndop = dop[P_LINDEX_ONE(i,j,k)] / vol[P_LINDEX_ONE(i,j,k)];
            else
              ndop = 0;
            
            alfa=ndop / (band.Ni * 2.0);
            
            if(alfa > 0.0)
              vadd_val[P_LINDEX_ONE(i,j,k)] = log( alfa + sqrt( 1.0 + alfa * alfa ));
            else if (alfa < 0)
              vadd_val[P_LINDEX_ONE(i,j,k)] = - log( sqrt( 1.0 + alfa * alfa ) - alfa );
            else vadd_val[P_LINDEX_ONE(i,j,k)] = 0;
          }

      if (p_jend == p_numy - 1)
        for (i = p_ibegin; i <= p_iend; i ++)
          for (k = p_kbegin; k <= p_kend; k ++)
                vadd_val[P_LINDEX_ONE(i,p_jend, k)] = vadd_val[P_LINDEX_ONE(i,p_jend - 1,k)];

      if (p_jbegin == 0)
        for (i = p_ibegin; i <= p_iend; i ++)
          for (k = p_kbegin; k <= p_kend; k ++)
                vadd_val[P_LINDEX_ONE(i,0, k)] = vadd_val[P_LINDEX_ONE(i,1,k)];


    for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
        for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
          charge_fac[P_LINDEX_ONE(i,j,k)] = vol[P_LINDEX_ONE(i,j,k)] * Nc / conc0;

    double Ez = ec0 * pow(hq0 * PI / tsi, 2) / (2 * mz[0] * em0 * spr0 * spr0);

    /* I assume that i = p_tox + 1 is in silicon, and vadd is the same every
    * where since the S/D contact is dopped equally everywhere*/
    if (p_jbegin == 0)  /* we are in charge of the source contact */
        fermi[0] = (-Vs ) / pot0;
    
    if (p_jend == p_numy - 1) /* we are in charge of drain contact */
        fermi[1] = (-Vd ) / pot0;
  
}

void MeshQuantities::init_epetra_map_vector() {

  int i, j,k;

  int *p_global_id1, *p_global_id2, *c_global_id1, *c_global_id2;
    
  /**
   * @brief
   * 参考了一下 3DMC manual 文档
   * 1. The X-direction is the horizontal direction from left to
   *  right along the length of the device.
   * 
   * 2. Y-direction is the direction along the length of the channel
   *  (FinFET 中沿着沟道长度的方向)
   * 3. Z-direction is the direction along the width of the channel.
   *  (FinFET 中平行与 Fin 高的方向)
   * 
   * 因此可以想象一下，以 FinFET 为例，沿着沟道长度的方向为Y轴
   *  这时，我们沿着 Y 轴将器件，按照 Y 轴上的点，得到一个个 XZ方向的切面
   */


  /**
   * @brief Epetra_Map class is used to define a mapping of global indices
   *  to local indices for distributed memory parallel computations.
   * 
   * @return a Pointer to a Epetra_Map object.
   * This pointer partitions the `c_numy` elements across the processors
   *  defined by the communicator `Comm`.
   * 
   * @param c_numy: The total number of cells in the Y-direction.
   * @param 0:      The global indices will start from 0.
   * @param *Comm:  A pointer to an Epetra_Comm object,
   *  which defines the communication context.
   * 
   * @example
   * 下面是一个例子，假设现在 Y方向上有100个 Cell，即 c_numy = 100
   *  运行时，使用4个 processors.
   *  此时，4 个 processors 中的信息分别为：
   * Processor 0 has 25 elements: 0 1 2 ... 24
   * Processor 1 has 25 elements: 25 26 27 ... 49
   * Processor 2 has 25 elements: 50 51 52 ... 74
   * Processor 3 has 25 elements: 75 76 77 ... 99
   * 此时，my_global_elements 存储的是当前 processor 中 cell 的 global_index
   * 而 num_local_elements 存储的是当前 processor 中的 cell 数目
   *  如，当 *Comm 代表 processor 0 时，my_global_elements 存储 [0~24]  
   *    而 num_local_elements 的值为 25。
   * 
   */
  cout << "epe_map_pos0" << endl;
  //cout << "c_numy: " << c_numy << endl;
  c_map_y = new Epetra_Map(c_numy, 0, *Comm);

  // Returns a pointer to an array containing the global indices 
  //  of the elements that are owned by the calling processor.

  cout << "epe_map_pos1" << endl;

  int * my_global_elements = c_map_y->MyGlobalElements();

  // Returns the number of elements that are owned by the calling processor.
  int num_local_elements = c_map_y->NumMyElements();

  cout << "epe_map_pos2" << endl;
    

  c_num_local_y = num_local_elements;     // 当前 processors 上，Y方向上的cell数目

  
  c_jbegin = my_global_elements[0];       // 当前 processor 上，y 方向上的起始 cell 的 global index
  c_jend = c_jbegin + c_num_local_y - 1;  // 当前 processor 上，y 方向上的最后一个 cell 的 global index

  p_jbegin = c_jbegin;  // 当前 Processor 上，y 方向上的起始 Point 的 global index
  p_jend = c_jend + 1;  // 当前 Processor 上，y 方向上的最后一个 Point 的 global index
  p_num_local_y = c_num_local_y + 1;  // 当前 Processor上，y 方向上的 Point 的 global index 数

  /**
   * 
   */
  p_num_local_nonoverlap_y = c_num_local_y; // 当前 processor上，y 方向上不重叠的Point 数目

  p_jbegin_nonoverlap = p_jbegin; // 当前 processor上，y 方向上不重叠的起始 Point 的 global index
  p_jend_nonoverlap = p_jend - 1; // 当前 processor上，y 方向上不重叠的最后一个 Point 的 global index
  

  /**
   * @details
   * p_jend 为 当前 Processor 上，Y 轴上的最后一个 Point 的 global index
   * p_numy 为整个 Device 中 Y 轴上的 Point 数
   * 当 global index == 整个 Device 上 Y 轴的 Point 数-1 时 
   *  说明此时只有一个 Processor 在运行，那么此时应该就 不需要 考虑 overlapping Cell
   * 
   * @details
   * 原先 Y 轴上的 Nonoverlap 的 Point = cell 数，也就是默认了最后一个 cell 为 Overlapping Cell
   *  但是，此时只有一个 Processor，不涉及通信，因此就无需考虑 Overlapping Cell
   *  那么，Y 轴上的所有 Point 都是 Nonoverlapping Point，
   *  那么，此时 Point 数 = Cell数 + 1
   * 
   * 由于此时最后一个 Point 也是 Nonoverlap 的 Point，那么最后一个 Point 的 gloabl index 也要重新加入进来
   * p_jend_nonoverlap 为 Y 轴上组成不重叠区域的最后一个 Point 的 Global index
   */
  if (p_jend == p_numy - 1) {   
    p_num_local_nonoverlap_y++; 
    p_jend_nonoverlap ++;
  }

  c_ibegin = 0;         // X 方向上的起始 Cell 的 Local Index
  c_iend = c_numx - 1;  // X 方向上的最后一个 Cell 的 Local Index
  c_isize = c_numx;     // X 方向上的 Cell 的 Local Index 数目
  p_ibegin = 0;         // X 方向上起始 Point 的 Local Index
  p_iend = c_numx;      // X 方向上最后一个 Point 的 Local Index
  p_isize = p_numx;     // X 方向上的 Point 的 Local Index 数目

  c_kbegin = 0;         // Z 方向上的起始 Cell 的 Local Index
  c_kend = c_numz - 1;  // Z 方向上的最后一个 Cell 的 Local Index
  c_ksize = c_numz;     // Z 方向上 Cell 的 Local Index 数目
  p_kbegin = 0;         // Z 方向上的起始 Point 的 Local Index
  p_kend = c_numz;      // Z 方向上的最后一个 Point 的 Local Index
  p_ksize = p_numz;     // Z 方向上的 Point 的 Local Index 数目

  c_numxz = c_numx * c_numz;    // XZ 面上的 Cell 数目
  p_numxz = p_numx * p_numz;    // XZ 面上的 Point 数目

  c_num_local = c_num_local_y * c_numxz;  // 整个 Device 中，属于当前 processor 的 cell 总数
  p_num_local = p_num_local_y  * p_numxz; // 整个 Device 中，属于当前 processor 的 Point 总数

  // 整个 Device 中，属于当前 Processor 的，不重叠的Region 的 Point 总数
  p_num_local_nonoverlap =p_num_local_nonoverlap_y * p_numxz;
    
  // 对当前 Processor，数  = 当前 Processor 上在Y方向上的 Cell 数
  // 这个数目是 cell数（including ghost cell) or just the number of ghost cell?
  //  我感觉更应该是 Total number of cell including ghost cell.
  c_num_local_y_ghost = c_num_local_y;
  c_jbegin_ghost = c_jbegin;    // 当前 Processor 上的起始 Ghost Cell 的 global index 
  c_jend_ghost = c_jend;        // 当前 Processor 上最后一个 Ghost Cell 的 global index
  

  /**
   * @details
   * [0] != 0，即第一个 cell 的 global index 不等于 0
   *  说明此时不是第一个 Processor
   * ! 那么，在这个 Processor 上，起始和最后一个cell 要各添加一个
   * 
   * @example
   * for example, 在 Processor 1 上，Cell分别为[3,4,5]
   * 此时，考虑 Ghost cell 时的 cell 为 [2,3,4,5,6] 
   * 这里的 2 就依靠 c_jbegin_ghost--
   * 而 6 就依靠 c_jend_ghost++
   */
  if (my_global_elements[0] != 0) {
    c_num_local_y_ghost ++;   // 此时，当前 Processor 上的 Ghost cell 数 + 1
    c_jbegin_ghost --;        // 此时，当前 Processor 上的 Ghost cell 的 global index 要-1
  }

  /**
   * @details
   * [num_local_elements - 1] != c_numy - 1 
   *  说明此时不是最后一个 Processor
   * ！那么，在这个 Processor 上，起始和最后一个 cell 都要各添加一个
   *  
   * * @example
   * for example, 在 Processor 1 上，Cell分别为[3,4,5]
   * 此时，考虑 Ghost cell 时的 cell 为 [2,3,4,5,6] 
   * 这里的 2 就依靠 c_jbegin_ghost--
   * 而 6 就依靠 c_jend_ghost++
   */
  if (my_global_elements[num_local_elements - 1] != c_numy - 1) {
    c_num_local_y_ghost ++; // y 方向上的 ghost 数
    c_jend_ghost ++;
  }
  
  // 当前 Processor 中第一个 ghost point 和 最后一个 ghost point
  p_jbegin_ghost = c_jbegin_ghost;
  p_jend_ghost = c_jend_ghost + 1;

#ifdef DEBUG
  cout << "c_jbegin = " << c_jbegin << endl
       << "c_jend = " << c_jend << endl
       << "c_jbegin_ghost = " << c_jbegin_ghost << endl
       << "c_jend_ghost = " << c_jend_ghost << endl;
#endif
    
  // 在考虑 Ghost cell 时，
  //  当前 Processor 上的 cell总数 和 point 总数
  // +1 是为了解决索引的问题
  // 例如，考虑 Ghost cell 时的 cell 为 [2,3,4,5,6] 
  // 那么 c_jend_ghost - c_jbegin_ghost = 6-2 = 4，而实际有5个 cell
  //  因此要 + 1
  c_num_local_ghost = (c_jend_ghost - c_jbegin_ghost + 1) * c_numxz;
  p_num_local_ghost = (c_jend_ghost - c_jbegin_ghost + 2) * p_numxz;

  cout << "epe_map_pos3" << endl;
  
  // 计算并存储每个 cell 的 global_index 
  c_global_id1 = (int *) malloc(sizeof(int) * c_num_local);
  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
      for (j = c_jbegin; j <= c_jend; j ++)
        c_global_id1[C_LINDEX_ONE(i,j,k)] = C_GINDEX(i,j,k);
  
  // 计算并存储 每个 ghost cell 的 global index
  c_global_id2 = (int *) malloc(sizeof(int) * c_num_local_ghost);
  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
      for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++)
	      c_global_id2[C_LINDEX_GHOST_ONE(i,j,k)] = C_GINDEX(i,j,k);
      
  cout << "epe_map_pos4" << endl;
  
  c_map = new Epetra_Map(-1, c_num_local, c_global_id1, 0, *Comm);

  free(c_global_id1);

  c_map_ghost = new Epetra_Map(-1, c_num_local_ghost, c_global_id2, 0, *Comm);

  free(c_global_id2);
  
  c_map_ghost_4 = new Epetra_Map(-1, c_num_local_ghost * 6, 0, *Comm);
  
  p_global_id1 = (int *) malloc(sizeof(int) * p_num_local_nonoverlap);
  
  p_global_id2 = (int *) malloc(sizeof(int) * p_num_local);

  cout << "epe_map_pos5" << endl;
  
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
	      p_global_id1[P_LINDEX_ONE(i,j,k)] = P_GINDEX(i,j,k);
    
  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin; j <= p_jend; j ++)
	      p_global_id2[P_LINDEX_ONE(i,j,k)] = P_GINDEX(i,j,k);

  cout << "epe_map_pos6" << endl;
  
  p_map_nonoverlap = new Epetra_Map(-1, p_num_local_nonoverlap, p_global_id1, 0, *Comm);
  
  p_map = new Epetra_Map(-1, p_num_local, p_global_id2, 0, *Comm);
  
  free(p_global_id1);
  free(p_global_id2);

  p_global_id1 = (int *) malloc(sizeof(int) * p_num_local_ghost);

  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_ghost; j <= p_jend_ghost; j ++)
	      p_global_id1[P_LINDEX_ONE_GHOST(i,j,k)] = P_GINDEX(i,j,k);
 
  p_map_ghost = new Epetra_Map(-1, p_num_local_ghost, p_global_id1, 0, *Comm);

  free(p_global_id1);

  cout << "epe_map_pos7" << endl;

  c_donor = new Epetra_Vector(*c_map_ghost);

  c_acceptor = new Epetra_Vector(*c_map_ghost);

  c_volume = new Epetra_Vector(*c_map_ghost);

  c_da = new Epetra_Vector(*c_map_ghost);
  
  c_material = new Epetra_IntVector(*c_map_ghost);

  c_desired_electron_number = new Epetra_IntVector(*c_map_ghost);
  
  c_desired_hole_number = new Epetra_IntVector(*c_map_ghost);
  
  c_InSurfRegion = new Epetra_IntVector(*c_map_ghost);

  c_quantumRegion = new Epetra_IntVector(*c_map_ghost);

  c_electron_num = new Epetra_IntVector(*c_map_ghost);

  c_hole_num = new Epetra_IntVector(*c_map_ghost);

  c_attached_contact = new Epetra_IntVector(*c_map_ghost);

  c_scatter_type = new Epetra_IntVector(*c_map_ghost);

  c_motion_rule = new Epetra_IntVector(*c_map_ghost_4);

  c_init_electron_charge = new Epetra_Vector(*c_map_ghost);

  c_init_hole_charge = new Epetra_Vector(*c_map_ghost);

  c_electron_charge = new Epetra_Vector(*c_map_ghost);

  c_hole_charge = new Epetra_Vector(*c_map_ghost);

  c_par_charge = new Epetra_Vector(*c_map_ghost);
  
  p_dop = new Epetra_Vector(*p_map);

  p_volume = new Epetra_Vector(*p_map);

  p_vadd = new Epetra_Vector(*p_map);
  
  p_par_charge = new Epetra_Vector(*p_map_nonoverlap);

  old_p_par_charge = new Epetra_Vector(*p_map_nonoverlap);
  
  par_list.resize(c_num_local_ghost);

  p_pot = new Epetra_Vector(*p_map);
//  p_qc_pot = new Epetra_Vector(*p_map);
  p_quantum_stat_pot = new Epetra_Vector(*p_map);

  p_poisson_pot = new Epetra_Vector(*p_map_nonoverlap);

  p_material = new Epetra_IntVector(*p_map_nonoverlap);

  p_poisson_pot_saved = new Epetra_Vector(*p_map_nonoverlap);
  
  new_Ec = new Epetra_Vector(*p_map_nonoverlap);

  p_nq = new Epetra_Vector(*p_map_nonoverlap);
  p_fermi_level = new Epetra_Vector(*p_map_nonoverlap);
  p_qc_fermi_level = new Epetra_Vector(*p_map_nonoverlap);

  p_rhs = new Epetra_Vector(*p_map_nonoverlap);
  /* used to fill the data from nonoverlap array to the overlap array */
  p_Importer = new Epetra_Import(*p_map, * p_map_nonoverlap);
  
  p_Identity_Importer = new Epetra_Import(*p_map_nonoverlap, *p_map_nonoverlap);

  //p_heat_Importer = new Epetra_Import(*p_map_ghost, *p_map_nonoverlap);

  /* used to fill the data in the ghost cells */
  c_Importer = new Epetra_Import(*c_map_ghost, *c_map);
  /* the following two cell variables do not include the ghost cells*/
  c_work1 = new Epetra_Vector(*c_map);

  c_work2 = new Epetra_Vector(*c_map);

  c_work3 = new Epetra_Vector(*c_map);

  p_work1 = new Epetra_Vector(*p_map_nonoverlap);

  c_field_x = new Epetra_Vector(*c_map_ghost);
  
  c_field_y = new Epetra_Vector(*c_map_ghost);

  c_field_z = new Epetra_Vector(*c_map_ghost);

  c_h_field_x = new Epetra_Vector(*c_map_ghost);
  
  c_h_field_y = new Epetra_Vector(*c_map_ghost);

  c_h_field_z = new Epetra_Vector(*c_map_ghost);

  stat_vxE = new Epetra_Vector(*p_map_nonoverlap);

  stat_vxH = new Epetra_Vector(*p_map_nonoverlap);

  stat_vyE = new Epetra_Vector(*p_map_nonoverlap);

  stat_vyH = new Epetra_Vector(*p_map_nonoverlap);

  stat_vzE = new Epetra_Vector(*p_map_nonoverlap);

  stat_vzH = new Epetra_Vector(*p_map_nonoverlap);

  stat_chargeE = new Epetra_Vector(*p_map_nonoverlap);

  stat_chargeH = new Epetra_Vector(*p_map_nonoverlap);

  stat_energyE = new Epetra_Vector(*p_map_nonoverlap);

  stat_energyH = new Epetra_Vector(*p_map_nonoverlap);

  stat_pot = new Epetra_Vector(*p_map_nonoverlap);

  stat_qc_pot = new Epetra_Vector(*p_map_nonoverlap);

  stat_pot_saved = new Epetra_Vector(*p_map_nonoverlap);

  stat_qc_pot_saved = new Epetra_Vector(*p_map_nonoverlap);

  stat_chargeE_saved = new Epetra_Vector(*p_map_nonoverlap);

  stat_chargeH_saved = new Epetra_Vector(*p_map_nonoverlap);

  stat_e_heat = new Epetra_Vector(*p_map_nonoverlap);

  stat_h_heat = new Epetra_Vector(*p_map_nonoverlap);

  stat_ec  = new Epetra_Vector(*c_map_ghost);
  
  stat_enum = new Epetra_Vector(*c_map_ghost);

  p_charge_fac = new Epetra_Vector(*p_map_nonoverlap);

  c_nearestSurf = new Epetra_IntVector(*c_map_ghost);
  c_EeffDirection = new Epetra_IntVector(*c_map_ghost);
  c_e_roughnessRate = new Epetra_Vector(*c_map_ghost);
  c_h_roughnessRate = new Epetra_Vector(*c_map_ghost);
  c_e_surfPhononRate = new Epetra_Vector(*c_map_ghost);
  c_h_surfPhononRate = new Epetra_Vector(*c_map_ghost);

  p_electron_heat = new Epetra_Vector(*p_map_ghost);

  p_hole_heat = new Epetra_Vector(*p_map_ghost);

  p_e_heat_weight = new Epetra_Vector(*p_map_ghost);

  p_h_heat_weight = new Epetra_Vector(*p_map_ghost);

  cout << "epe_map_pos8" << endl;
/*
  subbands.resize(p_num_local_nonoverlap_y);
  subband_vec.resize(p_num_local_nonoverlap_y);
  density_sub.resize(p_num_local_nonoverlap_y);

  for (i = p_jbegin_nonoverlap; i <= p_jend_nonoverlap; i ++) {
    subbands[P_SHIFT_NONOVERLAP_Y(i)].resize(valley_num);
    density_sub[P_SHIFT_NONOVERLAP_Y(i)].resize(valley_num);

    subband_vec[P_SHIFT_NONOVERLAP_Y(i)].resize(valley_num);
    for (j = 0; j < valley_num; j ++){
      subbands[P_SHIFT_NONOVERLAP_Y(i)][j].resize(max_subband);
      density_sub[P_SHIFT_NONOVERLAP_Y(i)][j].resize(max_subband);

      subband_vec[P_SHIFT_NONOVERLAP_Y(i)][j].resize(max_subband);
    }
  }
  */
}

void MeshQuantities::getInputData(char * FileName) {
  
  Trilinos_Util::InputFileReader fileReader(FileName);

  fileReader.ReadFile();
  
  grid_file_name = fileReader.Get("gridFile", "lgrid.txt");
  
  total_step = fileReader.Get("total_step", 100000);  //总步长

  /**
   * @brief The initial total number of Electrons in the whole device.
   * 
   * @attention This value MUST be set large enough 
   *  so that at least one electron exits in each grid cell according to electical neutrality. 
   *
   * Default Value: 100000 (1e6) 
   */
  electron_number = fileReader.Get("ElectronNumber", 100000);   
  
  /**
   * @brief The initial total number of Holes in the whole device.
   * 
   * @attention This value MUST be set large enough 
   *  so that at least one hole exits in each grid cell according to electical neutrality. 
   *
   * Default Value: 100000 (1e6) 
   */
  hole_number = fileReader.Get("HoleNumber", 100000);

  device_file_name = fileReader.Get("device_file_name", "ldg.txt");

  /**
  * @brief Step of printing the simulation info for debugging
  */
  debug_print_step= fileReader.Get("debug_print_step", 200);

  valley_num = fileReader.Get("valley", 3);   //波谷数

  max_subband = fileReader.Get("subband", 3);

  ml = fileReader.Get("mlong", 0.91);

  mt = fileReader.Get("mtrans", 0.19);

  mx[0] = mt; mx[1] = ml; mx[2] = mt;
  my[0] = mt; my[1] = mt; my[2] = ml;
  mz[0] = ml; mz[1] = mt; mz[2] = mt;

  Ncc = 2 * em0 * mt * eV0 / (PI * ec0 * hq0 * hq0);  //二维的Nc 实际上的公式 (2 * m* kB * T)/(pi * hbar^2)

  for (int i = 0;i < valley_num; i ++){
    Nccc[i] = 2 * em0 * sqrt(mx[i] * my[i]) * eV0 / (PI * ec0 * hq0 * hq0);
  }

  Flag_NonLinearPoisson = fileReader.Get("NonLinearPoisson", false);
                                     
  Flag_MultipleRefresh= fileReader.Get("Flag_MultipleRefresh", false);

  Rpar=2;
  Rsqr=2;
  pckill=1.0e-30;
  
    //RunControl
  Flag_QuantumCorrection = fileReader.Get("Flag_QuantumCorrection", false);
  
  Flag_LaterQC = fileReader.Get("Flag_LaterQC", false);
  
  Flag_SPE = fileReader.Get("Flag_SPE", false);

  quantum_start_step = fileReader.Get("quantum_start_step", 20);

  quantum_stat_step = fileReader.Get("quantum_stat_step", 20); 

  quantum_print_steps = fileReader.Get("quantum_print_steps", 1000); 

  /**
   * @brief Flag indicates if Surface Related Scatterings
   *  are considered in the simulation.
   * 
   * Default Value: true
   */
  Flag_SurfaceScatter = fileReader.Get("Flag_SurfaceScatter", true);
  
  /** 
   * @brief Flag indicates if Surface Roughness Scattering
   *  is considered in the simulation. 
   * 
   * Default Value: true
   */
  Flag_SurfaceRoughnessScatter = fileReader.Get("Flag_SurfaceRoughnessScatter", true);
  
  /**
   * @brief Flag indicates if Surface Phonon Scattering 
   *  is considered in the simulation.
   * 
   * Default Value: true
   */
  Flag_SurfacePhononScatter = fileReader.Get("Flag_SurfacePhononScatter", true);
  
  /**
   * @brief Flag indicates if Surface Impurity Scattering
   *  is considered in the simulation.
   * 
   * Default Value: true 
   */
  Flag_SurfaceImpurityScatter = fileReader.Get("Flag_SurfaceImpurityScatter", false);

  Flag_compute_potential = fileReader.Get("Flag_compute_potential", true);

  /**
   * @brief Flag indicates if the heat 
   *  is computed.
   * @attention Until now, I don't know whether the computation is self-consistently or not.
   * 
   * Default Value: false
   */
  Flag_compute_heat = fileReader.Get("Flag_compute_heat", false);
  
  Flag_calSurfscatt = fileReader.Get("Flag_calSurfscatt", false);
  
  Flag_test = fileReader.Get("Flag_test", false);

  /**
   * @brief Total number of time steps of considering heat generation
   * 
   * Default Value: 20000 (2e4) 
   */
  heat_steps = fileReader.Get("heat_steps", 20000);
  
  //  ndt = fileReader.Get("ndt", );
  
  dt = fileReader.Get("dt", 1e-16);

  Flag_restart = fileReader.Get("Flag_restart", false);

  restart_step = fileReader.Get("restart_step", 0);

  restart_filename = fileReader.Get("restart_filename", "nosuchname");

  stat_step = fileReader.Get("stat_step", 5000);

  stat_heat_step = fileReader.Get("stat_heat_step", 5000);
  
  qc_stat_step = fileReader.Get("qc_stat_step", 5);
  
  output_dir = fileReader.Get("OutPut_DIR", "../data");
  
  mr_step = fileReader.Get("mr_step", 10);

  char path_buf[256];
  
  memset(path_buf, 0, sizeof(path_buf));
  
  readlink("/proc/self/exe", path_buf, sizeof(path_buf));
  
  bs_path = dirname(path_buf) + string("/../input");

  kloem_table_file = fileReader.Get("kloem_file", "kloem.txt");
  kloab_table_file = fileReader.Get("kloab_file", "kloab.txt");
  ktoem_table_file = fileReader.Get("ktoem_file", "ktoem.txt");
  ktoab_table_file = fileReader.Get("ktoab_file", "ktoab.txt");
  klaem_table_file = fileReader.Get("klaem_file", "klaem.txt");
  klaab_table_file = fileReader.Get("klaab_file", "klaab.txt");

 
  
  //fileReader.ShowAll();
}

/** 
 * @brief 初始化物理参数，注意这个参数写死了，要换模型的话得自己改代码
 */
/**
 * @brief 初始化物理常量并做归一化，后续计算使用无量纲单位
 * @param filename 供读取温度等输入文件
 *
 * 流程：
 * 1) 设定材料开关（Si/GaAs）、质量密度、声速等本征参数；
 * 2) 读取器件温度，按 k_B T 构造能量、动量、长度、时间等基准尺度；
 * 3) 计算态密度 Nc、电流尺度 N_cur、量子势系数等归一化常数；
 * 4) 把介电常数、表面散射参数等转换到无量纲。
 */
void MeshQuantities::init_phpysical_parameter(char * filename) {
 
  frickel=1.0;          // Frickel 参数，保持默认 1
	
  sifl = true;          // 使用硅材料
  gaasfl = false;       // 未启用 GaAs

  psi_si = 4.05;        // Si 的功函数(eV)

  if(gaasfl)
    {
      //gaas bulk parameters
      //lattice constant
      sia0 = 5.64e-10;
      //mass density
      sirho= 5.36e3;
      //sound velocity
      siul =5.24e3;
      siut =2.47e3;
    }
  else if(sifl)
    {
      //si bulk parameters
      //lattice constant
      sia0 =5.43e-10;
      //mass density
      sirho=2.33e3;
      //sound velocity
      siul =9.05e3;
      siut =9.05e3;
    }
  //     get the material coefficients
  //____and normalize them
  read_device_input_temperature(filename); // 从输入文件读取温度
  T0 = device_temperature;
  cout << "Tem: " << T0 << endl;
  Tn=T0/300.0;
  // 初始参数（后续可改为从 input 文件读取）
  double ml_val = 0.916;   // 纵向有效质量（相对 m0）
  double mt_val = 0.19;    // 横向有效质量（相对 m0）
  double alpha_val = 0.5;  // 非抛物线性系数 (1/eV)
  double alpha;
  double alpha_norm;
	
  //     energy [eV]/electon rest mass [kg]/Planck's constant [eVs] /
  //____electron charge [As]
  eV0   =BOLTZ*T0;   // 能量尺度 k_B T [eV]
  em0   =EM;         // 电子质量 [kg]
  hq0   =PLANCK;     // 普朗克常量 [eV*s]
  ec0   =EC;         // 电子电荷 [C]
  //     momentum [eVs/m]/r-space [m]/k-space [1/m]/time [s] /
  //____velocity [m/s]
  rmom0 =sqrt((em0/ec0)*eV0);  // 动量尺度
  spr0  =hq0/rmom0;           // 长度尺度
  spk0  =1.0/spr0;            // 波矢尺度
  time0 =hq0/eV0;            // 时间尺度 (h/k_BT)
  velo0 =spr0/time0;         // 速度尺度
  cvr   =CLIGHT/velo0;       // 光速对应的无量纲数
  //____el. potential [V]/el. field [V/m]/concentration [1/m**3]
  pot0  =eV0;                // 电势尺度
  field0=pot0/spr0;          // 电场尺度
  conc0 =spk0*spk0*spk0;     // 浓度尺度
  //____mass density [kg/m**3]
  dens0 =em0*conc0;          // 质量密度尺度

  // 归一化有效质量与非抛物线参数
  ml = ml_val;
  mt = mt_val;
  alpha = alpha_val * eV0;
  alpha_norm = alpha;

  //____deformation potential constant [eV/m]/scattering rate [1/s]
  dpc0  =field0;            // 应变势尺度
  scrt0 =1.0/time0;         // 散射率尺度 k_BT/h  

  //____current [A/m]
  curr0 =ec0/time0;

  //effective density of state: 2.82 x 1e25 , aproximately
  Nc = 2 * pow(1.08 * em0 * pot0 / (hq0 * hq0 * ec0 * 2 * PI), 1.5); // 有效态密度

  N_cur = 2 * hq0 / (em0 * PI * PI) * pow((em0 * eV0) / (hq0 * hq0), 2) / ec0; // 电流尺度 A/m^2

  //quantum potential coefficient
  QuantumPotentialCoef=hq0 * hq0 * ec0 / (12.0 * 0.26 * em0 * spr0 * spr0 * pot0); // 量子势系数
     

  //____Si bulk parameters
  sia0=sia0/spr0;      // 归一化晶格常数
  sirho=sirho/dens0;   // 归一化质量密度
  siul=siul/velo0;     // 归一化纵向声速
  siut=siut/velo0;     // 归一化横向声速

  //____temperature dependent band gap
  if(T0<190.0)
    sieg=(1.170+1.059e-5*T0-6.05e-7*T0*T0)/eV0;
  else if(T0<250.0)
    sieg=(1.17850-9.025e-5*T0-3.05e-7*T0*T0)/eV0;
  else
    sieg=(1.2060-2.730e-4*T0)/eV0;

  //____Defines a0pi
  a0pi=TWOPI/sia0;

  //____get relative dielectric constant and normalize (eps_vacuum <> 1 within
  //     the program)
  eps[VACUUM]= 1.00/(4*PI*cvr*FSC);
  eps[OXIDE]= 3.90/(4*PI*cvr*FSC);
  if(gaasfl)
    eps[SILICON]=12.90/(4*PI*cvr*FSC);
  else if (sifl)
    eps[SILICON]=11.70/(4*PI*cvr*FSC);
  SurfSc_ail=2.0e-9/spr0;//m
  SurfSc_delta=0.1126755e-9/spr0;//m
  SurfSc_XIph=12.0;     //eV
  SurfSc_theta=2.857703;
  SurfSc_gama=1.5;
  SurfSc_Nbmod=2.0;
  SurfSc_Npf=-0.6;
  SurfSc_Pft=0.7;
  SurfSc_Pfn=0.3;
  SurfSc_Kpha=3600e-4;//A*s*s/kg
  SurfSc_Kt=-2.24;
  SurfSc_Fsimp=0.6;
  SurfSc_Rshmin=0.5;
  SurfSc_Rshmax=10;
  SurfSc_GAMMAn=0.08;
}

/*scale the input: dopping concentration, potential, length*/
void MeshQuantities::scaling() {

  int i, j, k;
  double * donor, * acceptor;
  
  for (i = 0;i < dx.size(); i ++)
    dx[i] /= spr0;
  for (i = 0;i < dy.size(); i ++)
    dy[i] /= spr0;
  for (i = 0;i < dz.size(); i ++)
    dz[i] /= spr0;

  for (i = 0;i < lx.size(); i ++)
    lx[i] /= spr0;
  for (i = 0;i < ly.size(); i ++)
    ly[i] /= spr0;
  for (i = 0;i < lz.size(); i ++)
    lz[i] /= spr0;

  dt /= time0;
  
  inject_const = N_cur * dt * time0 * spr0 * spr0; 

  tsi = tsi * 1e-9 / spr0;

  for (i = 0;i < contact.size(); i ++){
      contact[i].CurrentVapp /= pot0;
    contact[i].PhiMS /= pot0;
  }

  c_donor->ExtractView(&donor);

  c_acceptor->ExtractView(&acceptor);
  
  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
      for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++)
      {
        donor[C_LINDEX_GHOST_ONE(i,j,k)] /= conc0;
        acceptor[C_LINDEX_GHOST_ONE(i,j,k)] /= conc0;
      }
}

int MeshQuantities::distance2index(double pos, vector<double> &cut) {
  
  for (int i = 0;i < cut.size(); i ++)
    if (fabs(pos - cut[i] * 1e9) < 1e-4) 
      return i;

  cout << pos << " not aligned " << endl;
  exit(1);
}

void MeshQuantities::get_cube_range(ifstream &ifile, double * pos, int * index) {
  int i;
  for (i = 0;i < 6; i ++)
    ifile >> pos[i] ;
  index[0] = distance2index(pos[0], lx);
  index[1] = distance2index(pos[1], lx);
  index[2] = distance2index(pos[2], ly);
  index[3] = distance2index(pos[3], ly);
  index[4] = distance2index(pos[4], lz);
  index[5] = distance2index(pos[5], lz);
}

void MeshQuantities::read_device_file() {

  cout << " read device file" << endl;
  
  int index[6];
  double cube_pos[6];
  string dir_name, motion_name[6], region_name;
  map<string, int> dir, motion, op, region_type;
  string para_type;
  
  // for The motionplane Command
  // The direction in which the particle hits the plane.
  dir["UP"] = UP;         // The negative X direction
  dir["DOWN"] = DOWN;     // The positive X direction
  dir["LEFT"] = LEFT;     // The negative Y direction
  dir["RIGHT"] = RIGHT;   // The positive Y direction
  dir["FRONT"] = FRONT;   // The negative Z direction
  dir["BACK"] = BACK;     // The positive Z direction

  /* ---------------------------------------------------------------------- */
  /**
   * Motion rule at the upper boundaries of the cells within the region.
   */
  /**
   * @brief Particles step through the interface and resume propagation.
   */
  motion["PASS"] = PASS;

  /**
   * @brief The particle is reflected at the interface and does not change the cell. 
   *  The velocity component perpendicular to the interface is reversed.
   */
  motion["REFLECT"] = REFLECT;

  /**
   * @brief The particle hits the Si/SiO2 interface.
   *  The particle is either reflected or diffusively scattered and does not change the cell. 
   */
  motion["SCATTOX"] = SCATTOX;

  /**
   * @brief A copy of the original particle is generated and placed on the other side of the interface. 
   *  The the PERIOD boundary conditions is applied to the original particle.
   */
  motion["GENERATE"] = GENERATE;

  /** 
   * @brief The particle is destoryed. Used for contacts. 
   */
  motion["CATCH"] = CATCH;

  /**
   * @brief The particle is moved to the opposite side of the surface of the cell.
   *  Only the component of the real space vector which is perpendicular to the surface is changed.
   */
  motion["PERIOD"] = PERIOD;

  /**
   * @brief Similar to GENERATE, except that the original particle is applied with REFLECT boundary condition. 
   */
  motion["GENREF"] = GENREF;

  /**
   * I'm not sure.
   * @brief The particle is destoryed by the Gate?
   */
  motion["CATCHGATE"] = CATCHGATE;

  double tmp, tmp1, tmp2, tmp3, tmp4, tmp5, delta;
  int itmp, ktmp, inum, knum, eig_num;
  int surfnum, surftype, surfdir;
  double surfpos;
  

  // op: different Commands
  /**
   * @brief 指定 cuboid region 为 n-type 掺杂。
   *  仅为 Silicon Region 设置。
   *  
   * 当区域重叠时，最后一个命令将会覆盖前一个命令
   * 
   * @details 常见的使用方法：donor xbegin xend ybegin yend zbegin zend conc
   *  例：donor 0 40  -32 -10  0 10  1e20 
   */
  op["donor"] = 1;

  /**
   * @brief 指定 cuboid region 为 p-type 掺杂。
   *  仅为 Silicon Region 设置。
   *  
   * 当区域重叠时，最后一个命令将会覆盖前一个命令
   * 
   * @details 常见的使用方法：acceptor xbegin xend ybegin yend zbegin zend conc
   *  例：acceptor 0 40 -10 10 0 10 1e18 
   */   
  op["acceptor"] = 2;

  /**
   * @brief 此命令允许指定网格中指定长方体区域的材质。
   * 如果网格的一个区域没有被器件占据，它必须被明确地指定为真空。电场和粒子都不能进入这些区域。
   * 这些区域的边界与网格线对齐。
   * 可以使用许多区域语句，在出现重叠区域的情况下，最后一个语句将覆盖前面的区域语句。
   * 
   * @details 使用方法：region xbegin xend ybegin yend zbegin zend mat_type
   *  例： region  0 60 -32 32 0 10 SILICON
   */
  op["region"] = 3;

  /**
   * @brief 每次粒子到达一个 Cell 的边界时，粒子的行为被打断，需要设置合适的行为。
   *  Cell 的边界为 网格的边界线。
   * 
   * @details 使用方法： motioncube xbegin xend ybegin yend zbegin zend rul_up rul_down rul_left rul_right rul_front rul_back
   * 例： motioncube   -130  39 -46  46 -19  19     PASS      PASS      PASS      PASS      PASS      PASS
   * 在上面的例子中，后面的6个pass代表在 X,Y,Z 的正负（共6个）方向上遇到边界时的运动方式。PASS 表示粒子会穿过边界，继续运动。
   */
  op["motioncube"] = 4;

  /**
   * @brief attachcontact 命令用于指定硅中的哪个区域属于 Contact
   *  
   * @attention 该区域必须包含平行与 Contact 平面的2层 Cell。
   *  此外，这个命令来设置并不能用于分配 gate Contact。
   *  
   * @details
   *  1. 每个 Contact 都对应一行 attachcontact 命令，并需要提供一个编号（大于0），不同的接触有不同的索引值。
   *  2. 由于程序中没有初始化 Silicon Contact 的代码，
   *    因此在设置器件的 Contact 前，用户必须将器件中的全部 silicon 区域的索引值设置为 0
   * 
   * @details attachcontact xbegin xend ybegin yend zbegin zend index
   * 例： attachcontact -130  39 -46  46 -19  19        0
   *     attachcontact   35  39 -46 -16 -14  14        1
   *  其中，第一行代表将整个器件中的 silicon 区域的索引值设置为 0 
   *       第二行才是设置实际的 contact
   */
  op["attachcontact"] = 5;

  /**
   * @brief contact 命令用于为器件施加电压
   *  contact 都是二维的平面，只能被垂直放置于 X-, Y- 和 Z- 方向。
   *  由 attachcontact 分配的 Silicon 中的 Contact 以及器件自身的 Gate Contact 都能用 Contact 来设置电压。
   * 
   * @details
   *  contact
   *  N phims (这是2个值)                            N代表contact的平面数目；phims为Workfunction difference between semiconductor and metal. 主要用于 gate contact
   *  xbegin1 xend1 ybegin1 yend1 zbegin1 zend1    1~N 不同的 Boundary
   *  xbegin2 xend2 ybegin2 yend2 zbegin2 zend2
   *  ...
   *  xbeginN xendN ybeginN yendN zbeginN zendN
   *  voltage                                      施加的电压值
   * 
   * 两个例子：
   * (1)
   * contact
   * 1 0
   * 39  39 -46 -16 -14  14
   * 0.0
   * 
   * (2)
   * contact
   * 3 -0.33
   * 0  35 -10  10  -5  -5
   * 0  35 -10  10   5   5
   * 35  35 -10  10  -5   5
   * 1.0
   * 
   * @details 便于理解的话，可以想象，S/D 这类的接触，都是薄薄一层，所以只有一个边界
   * 而 Gata Contact 这类接触，都要考虑 SiO2 的厚度，因此需要考虑不同方向上的边界。
   */
  op["contact"] = 6;

  /**
   * @brief motionplane 命令用于设置 粒子在某个平面的特定方向上的运动规则
   * motionplanes 命令仅能用于垂直于 X-, Y-, Z- 方向的平面
   * 
   * @details motionplane xbegin xend ybegin yend zbegin zend dir mot_rule
   * 其中 dir 代表粒子接触平面的方向，mot_rule 代表运动的规则
   * 
   * 来看两个例子：
   * (1) motionplane    37  37 -46 -16 -14  14     DOWN     CATCH
   * (2) motionplane    37  37 -46 -16 -14  14       UP    GENREF
   * 
   * 根据程序中设置 X，Y，Z方向的规则，可以看出 37  37 -46 -16 -14  14 是一个 YZ 的平面
   * 因此，两个例子分别给出了，当粒子碰撞到这个平面的上下表面时的运动规则。
   */
  op["motionplane"] = 7;

  /**
   * @brief 
   * 反正有2个参数，
   * 一个为 c_desired_electron_number （像是在 单个 cell 中理想电子数目）
   * 另一个为 c_desired_hole_number   （像是在 单个 cell 中理想的空穴数目）
   */
  op["parnumber"] = 8;

  /**
   * @brief 为 Multiple Refresh 方法提供粒子数
   * 
   * @details 使用方法：default_par_number ele_num hole_num
   * 
   */ 
  op["default_par_number"] = 9;

  // 只有一个参数，参数为 c_scatter_type, 是否代表了 单个 cell 中的散射类型？
  op["ScatterArea"] = 12;

/* VsVdVg and Vgrange are not used anymore */
  op["VsVdVg"] = 13;  // 不再使用
  op["vgrange"] = 14; // 不再使用

  /**
   * 
   * @details 如果用户不考虑使用自己的参数，则应该使用默认值：ep_parm 4 0.5 4 0.5 4 0.5 3.1
   */
  op["ep_parm"] = 15;

  /**
   * @brief 设置 Si/SiO2 界面的相关信息
   * 
   * @details 使用方法如下：
   * surfaces N                       // N 代表 Si/Oxide 界面的数目        
   * surf_type1 surf_pos1 surf_dir1   // type: 界面的类型，0 代表垂直于X轴，1 代表垂直于Y轴，2 代表垂直于Z轴
   * surf_type2 surf_pos2 surf_dir2   // pos:  代表界面的坐标
   * …
   * surf_typeN surf_posN surf_dirN   // dir:  界面的方向。表面的方向被定义为硅中的粒子向氧化物垂直移动的方向。共有6个方向
   * 
   * 
   */
  op["surfaces"] = 16;

  /**
   * @brief 该命令指定了发生 surface scattering 的 region，包括 surface roughness 和 surface phonon
   * 
   * @details surface_scatter_range xbegin xend ybegin yend zbegin zend
   */
  op["surface_scatter_range"] = 17;

  /**
   * @brief 该命令指定了考虑 effective potential quantum correction 的区域
   * 
   * @details quantumRegion xbegin xend ybegin yend zbegin zend
   */
  op["quantumRegion"] = 18;
  
  region_type["VACUUM"] = VACUUM;
  region_type["OXIDE"] = OXIDE;
  region_type["SILICON"] = SILICON;
  region_type["IGZO"] = IGZO;

  ifstream ifile;
  
  ifile.open(device_file_name.c_str());
  ifile >> para_type;
  
  cout << para_type << endl;

  user_cmd cmd;

  while (para_type != "end") {
    cmd.type = op[para_type];

    switch (op[para_type]){
    case 1:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> tmp;
      cmd.dbl_param[0] = tmp;
      break;
    case 2:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> tmp;
      cmd.dbl_param[0] = tmp;
      break;
    case 3:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> region_name;
      cmd.int_param[0] = region_type[region_name];
      if (cmd.int_param[0] == SILICON){
        cmd.int_param[1] = default_electron_num;
        cmd.int_param[2] = default_hole_num;
      }
      else if (cmd.int_param[0] == IGZO){
        cmd.int_param[1] = default_electron_num;
        cmd.int_param[2] = default_electron_num;
      } 
      else {
        cmd.int_param[1] = cmd.int_param[2] = 0;
      }
      break;
    case 7:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> dir_name >> motion_name[0] ;
      cmd.int_param[0] = dir[dir_name];
      cmd.int_param[1] = motion[motion_name[0]];
      break;
    case 4:
      get_cube_range(ifile, cube_pos, cmd.range);
      for (int i = 0;i < 6; i ++){
        ifile >> motion_name[i];
        cmd.int_param[i] = motion[motion_name[i]];
      }
      break;
    case 5:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> itmp;
      cmd.int_param[0] = itmp;
      break;
    case 8:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> itmp ;
      cmd.int_param[0] = itmp;
      ifile >> itmp;
      cmd.int_param[1] = itmp;
      break;
    case 9:
      ifile >> default_electron_num;
      ifile >> default_hole_num;
      break;
    case 6:
      Contact cont;
      cont.reset();
      ifile >> cont.NumContactPlane
          >> cont.PhiMS;
      for (int i = 0;i < cont.NumContactPlane; i ++){
        get_cube_range(ifile, cube_pos, index);
        cont.BeginI[i] = index[0];
        cont.EndI[i] = index[1];
        cont.BeginJ[i] = index[2];
        cont.EndJ[i] = index[3];
	      cont.BeginK[i] = index[4];
        cont.EndK[i] = index[5];
      }
      ifile >> cont.CurrentVapp;

      contact.push_back(cont);
      
      break;

    case 12:
      get_cube_range(ifile, cube_pos, cmd.range);
      ifile >> itmp;
      cmd.int_param[0] = itmp;
      break;
    case 13:
      ifile >> Vs >> Vd >> Vg >> phi_top;
      Eg = -Vg + phi_top - psi_si;
      break;
    case 14:
      ifile >> tmp >> tmp1 >> tmp2 >> tmp3;
      p_gbegin = distance2index(tmp, ly);
      p_gend = distance2index(tmp1, ly);
      p_tox = distance2index(tmp2, lx);
      p_box = distance2index(tmp3, lx);
      tsi = tmp3 - tmp2;

      if ((tsi < 2) && (T0 > 250))
	        fermi_order = 0;
      else 
          fermi_order = 0.5;
      break;
    case 15:
      ifile >> qc_xratio >> qc_xtheta >> qc_yratio >> qc_ytheta >> qc_zratio >> qc_ztheta >> Eb;

      qc_xtheta /= 1e9 * spr0; 
      qc_ytheta /= 1e9 * spr0; 
      qc_ztheta /= 1e9 * spr0; 
      Eb /= pot0;
      break;
    case 16:  // Keyword: surfaces
      ifile >> NumSurface;  // Total number of Si/Oxide surfaces in the devices.
      for (int i = 0;i < NumSurface; i ++){
        ifile >> surftype >> surfpos >> surfdir;
        SurfaceType[i] = surftype;  // 0: yz surface; 1: xz surface; 2: xy surface
        SurfacePosition[i] = surfpos * 1e-9 / spr0;  // 位置，转成内部无量纲长度
        SurfaceDir[i] = surfdir;  // 法向方向（粒子进氧化层的方向）
      }
      break;
    case 17:
      get_cube_range(ifile, cube_pos, cmd.range);
      break;
    case 18:
      get_cube_range(ifile, cube_pos, cmd.range);
      break;
    // default:
    //   if (mpi_rank == 0)
    //     cout << "unrecognized option: " << para_type << endl;
    }
    cmd_list.push_back(cmd);
    ifile >> para_type;
  }

  /* this is not used anymore */
/*  if (mpi_rank == 0)
    cout  << "Vs = "  << Vs << ',' << "Vd = " << Vd << ',' << "Vg = " << Vg << endl;
    */
}

bool file_exists(const std::string& filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
}

/* load the mesh file */
/**
 * @brief 读取 lgrid.txt 文件（存储了网格信息）中的信息
 */
void MeshQuantities::read_grid_file() {
  
  ifstream ifile;
  double pos;


    // Print the current working directory for debugging purposes
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        std::cerr << "getcwd() error: " << strerror(errno) << std::endl;
    }

    // Construct the full path to the grid file (prefer current working dir)
    std::string full_path = grid_file_name;
    if (!file_exists(full_path)) {
      // try relative to finfet/ because historical paths used that layout
      full_path = std::string("finfet/") + grid_file_name;
    }
    std::cout << "Attempting to open file: " << full_path << std::endl;

    // Check if the file exists
    // Why?????????????
    if (!file_exists(full_path)) {
      std::cerr << "File does not exist: " << full_path << std::endl;
      //return; // Exit the function if the file does not exist
    }

    // Check file permissions
    struct stat file_stat;
    if (stat(full_path.c_str(), &file_stat) != 0) {
      std::cerr << "Failed to get file status: " << strerror(errno) << std::endl;
      return;
    }

    if (!(file_stat.st_mode & S_IRUSR)) {
      std::cerr << "File does not have read permissions: " << full_path << std::endl;
      return;
    }

    

    

  ifile.open(full_path.c_str());

    if (ifile) {
          std::cout << "File opened successfully: " << full_path << std::endl;
      } else {
          std::cerr << "Failed to open file: " << full_path << std::endl;
          std::cerr << "Error: " << strerror(errno) << std::endl;
          return; // Exit the function if the file cannot be opened
      }
  
  // p_numx: the number of coordinates in X- direction.
  // 注意：它代表了X方向上一行的点的个数。
  ifile >> p_numx;

  // each line consists a coordinate.
  for (int i = 0 ; i < p_numx; ++ i) {
    ifile >> pos; // X 方向上的每个坐标
    lx.push_back(pos * 1e-6);
    if (i > 0) dx.push_back(lx[i] - lx[i - 1]);
  }
  
  // p_numy: the number of coordinates in Y- direction.
  ifile >> p_numy;
  
  // each line consists a coordinate.
  for (int i = 0 ; i < p_numy; ++ i){
    ifile >> pos; // Y 方向上的每个坐标
    ly.push_back(pos * 1e-6);
    if (i > 0) dy.push_back(ly[i] - ly[i - 1]);
  }
  
  // p_numz: the number of coordinates in Z- direction.
  ifile >> p_numz;

  for (int i = 0 ; i < p_numz; ++ i){
    ifile >> pos; // Z 方向上的每个坐标
    lz.push_back(pos * 1e-6);
    if (i > 0) dz.push_back(lz[i] - lz[i - 1]);
  }
 
  ifile.close();

  // The number of Cells in the X-, Y- and Z- direction.
  c_numx = p_numx - 1;
  c_numy = p_numy - 1;
  c_numz = p_numz - 1;

  std::cout << "c_numx: " << c_numx << "  "
          << "p_numx: " << p_numx << "  "
          << "c_numy: " << c_numy << "  "
          << "p_numy: " << p_numy << "   "
          << "c_numz: " << c_numz << "  "
          << "p_numz: " << p_numz << endl;
}

Epetra_CrsMatrix * MeshQuantities::init_matrix(){
   int l = 0;
   Epetra_CrsMatrix * mat;
   int indices[7];
   double values[7];
   int NumEntries;
   int i,j,k;
   double eps0;
   double AA, BB, CC, DD;
   int flag, icont, iplane;
   int *material;

   vector<vector<double> > nonzero;
   vector<int> entry_num;

   nonzero.resize(p_num_local_nonoverlap);

   for (i = 0;i < p_num_local_nonoverlap; i ++)
     for (j = 0;j < 7; j ++)
       nonzero[i].push_back(0);
   
   c_material->ExtractView(&material);

   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin_ghost; j <= c_jend; j ++)
       if (material[C_LINDEX_GHOST_ONE(i,j,k)] != VACUUM) {

	 eps0 = eps[material[C_LINDEX_GHOST_ONE(i,j,k)]];
	 AA = -0.25 * eps0 * dx[i] * dy[j] / dz[k];
	 BB = -0.25 * eps0 * dx[i] * dz[k] / dy[j];
	 CC = -0.25 * eps0 * dz[k] * dy[j] / dx[i];
	 DD = -AA - BB - CC;
          
	 if (j >= p_jbegin_nonoverlap){
	   nonzero[P_LINDEX_ONE(i,j,k)][5] += AA;
	   nonzero[P_LINDEX_ONE(i,j,k)][3] += BB;
	   nonzero[P_LINDEX_ONE(i,j,k)][1] += CC;
	   nonzero[P_LINDEX_ONE(i,j,k)][6] += DD;

	   nonzero[P_LINDEX_ONE(i,j,k + 1)][4] += AA;
	   nonzero[P_LINDEX_ONE(i,j,k + 1)][3] += BB;
	   nonzero[P_LINDEX_ONE(i,j,k + 1)][1] += CC;
	   nonzero[P_LINDEX_ONE(i,j,k + 1)][6] += DD;

	   nonzero[P_LINDEX_ONE(i + 1,j,k)][5] += AA;
	   nonzero[P_LINDEX_ONE(i + 1,j,k)][3] += BB;
	   nonzero[P_LINDEX_ONE(i + 1,j,k)][0] += CC;
	   nonzero[P_LINDEX_ONE(i + 1,j,k)][6] += DD;

	   nonzero[P_LINDEX_ONE(i + 1,j,k + 1)][4] += AA;
	   nonzero[P_LINDEX_ONE(i + 1,j,k + 1)][3] += BB;
	   nonzero[P_LINDEX_ONE(i + 1,j,k + 1)][0] += CC;
	   nonzero[P_LINDEX_ONE(i + 1,j,k + 1)][6] += DD;
	 }

	 if (j + 1 <= p_jend_nonoverlap){
	   nonzero[P_LINDEX_ONE(i,j + 1,k)][5] += AA;
	   nonzero[P_LINDEX_ONE(i,j + 1,k)][2] += BB;
	   nonzero[P_LINDEX_ONE(i,j + 1,k)][1] += CC;
	   nonzero[P_LINDEX_ONE(i,j + 1,k)][6] += DD;

	   nonzero[P_LINDEX_ONE(i,j + 1,k + 1)][4] += AA;
	   nonzero[P_LINDEX_ONE(i,j + 1,k + 1)][2] += BB;
	   nonzero[P_LINDEX_ONE(i,j + 1,k + 1)][1] += CC;
	   nonzero[P_LINDEX_ONE(i,j + 1,k + 1)][6] += DD;

	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k)][5] += AA;
	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k)][2] += BB;
	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k)][0] += CC;
	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k)][6] += DD;

	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k + 1)][4] += AA;
	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k + 1)][2] += BB;
	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k + 1)][0] += CC;
	   nonzero[P_LINDEX_ONE(i + 1,j + 1,k + 1)][6] += DD;
	 }
       }

   for (i = 0;i < p_num_local_nonoverlap; i ++){
     flag = false;
     for (j = 0;j < 7; j ++)
       if (nonzero[i][j] != 0) {
	 flag = true;
	 break;
       }
     if (!flag) nonzero[i][6] = 1;
   }
 
   for(icont=0;icont<contact.size();icont++)
    for(iplane=0;iplane<contact[icont].NumContactPlane;iplane++)
      for(i=contact[icont].BeginI[iplane];i<=contact[icont].EndI[iplane];i++)
	for(k=contact[icont].BeginK[iplane];k<=contact[icont].EndK[iplane];k++)
	  for(j=contact[icont].BeginJ[iplane];j<=contact[icont].EndJ[iplane];j++)
	    if ((j >= p_jbegin_nonoverlap) && (j <= p_jend_nonoverlap))
            {
	      for (l = 0;l < 7; l ++)
		nonzero[P_LINDEX_ONE(i,j,k)][l] = 0;
	      nonzero[P_LINDEX_ONE(i,j,k)][6] = 1;

            }
/*
    for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      if (!NOT_GHOST_BC(j)) 
	for (i = p_ibegin; i <= p_iend; i ++)
	  for (k = p_kbegin; k <= p_kend; k ++){
	    for (l = 0;l < 7; l ++)
	      nonzero[P_LINDEX_ONE(i,j,k)][l] = 0;
	    nonzero[P_LINDEX_ONE(i,j,k)][6] = 1;
      }
      */


   entry_num.resize(p_num_local_nonoverlap);

   for (i = 0;i < p_num_local_nonoverlap; i ++){
     entry_num[i] = 0;
     for (j = 0;j < 7; j ++)
       if (nonzero[i][j] != 0) 
	 entry_num[i] ++;
   }	

  mat = new Epetra_CrsMatrix(Copy, *p_map_nonoverlap, &(entry_num[0]));

  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++){

	NumEntries = 0;

	if (nonzero[P_LINDEX_ONE(i,j,k)][0] != 0){
	  indices[NumEntries] = P_GINDEX(i - 1,j,k);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][0]; 
	  NumEntries ++;
	}
	if (nonzero[P_LINDEX_ONE(i,j,k)][1] != 0){
	  indices[NumEntries] = P_GINDEX(i + 1,j,k);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][1]; 
	  NumEntries ++;
	}
	if (nonzero[P_LINDEX_ONE(i,j,k)][2] != 0){
	  indices[NumEntries] = P_GINDEX(i ,j - 1,k);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][2]; 
	  NumEntries ++;
	}
	if (nonzero[P_LINDEX_ONE(i,j,k)][3] != 0){
	  indices[NumEntries] = P_GINDEX(i ,j + 1,k);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][3]; 
	  NumEntries ++;
	}
	if (nonzero[P_LINDEX_ONE(i,j,k)][4] != 0){
	  indices[NumEntries] = P_GINDEX(i,j,k - 1);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][4]; 
	  NumEntries ++;
	}
	if (nonzero[P_LINDEX_ONE(i,j,k)][5] != 0){
	  indices[NumEntries] = P_GINDEX(i ,j,k + 1);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][5]; 
	  NumEntries ++;
	}
	if (nonzero[P_LINDEX_ONE(i,j,k)][6] != 0){
	  indices[NumEntries] = P_GINDEX(i,j,k);
	  values[NumEntries] = nonzero[P_LINDEX_ONE(i,j,k)][6]; 
	  NumEntries ++;
	}
     mat->InsertGlobalValues(P_GINDEX(i,j,k), NumEntries, values, indices);

    }

  //cout << *mat ;
   mat->FillComplete();
   /*while debugging, you can output the matrix element */

   return mat;
}

/* ----------------------------------------------------------------------- */
/**
 * @brief 
 *  Add values from neighboring process to the ghost cells of the local process,
 *    ensuring that the ghost cells contain the correct accumulated values from 
 *    the adjacent process.
 * 
 * @param p_vec         Pointer to the Epetra_Vector that holds the values to be reduced
 * @param ghost_width   The width of the ghost cells.
 */
void MeshQuantities::reduce_add(Epetra_Vector * p_vec, int ghost_width) {

  // Used for non-blocking MPI Communication
  MPI_Request request;
  MPI_Status status;

  double * p_vec_val,   // Pointer to the values 
         * send_buf,    // Buffers for sending data
         * recv_buf;    // Buffers for receiving data

  int send_jnum, send_num, recv_num, recv_jnum;

  int jbegin_index, jend_index, i,j,k;

  /* --- Determine Range of Cells --- */
  if (mpi_rank > 0)
    jbegin_index = p_jbegin - ghost_width;
  else 
    jbegin_index = p_jbegin;

  if (mpi_rank < mpi_size - 1)
    jend_index = p_jend + ghost_width;
  else 
    jend_index = p_jend;

  // Extract the view of the vector values into p_vec_val
  p_vec->ExtractView(&p_vec_val);


  /* --- Communication with Neighboring Processes --- */
  // Sending data to the next mpi_rank
  if (mpi_rank != mpi_size - 1) 
  {
    send_jnum = jend_index - p_jend + 1;

    send_num = send_jnum * p_numxz;

    send_buf = (double *) malloc(sizeof(double) * send_num);
  
    for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
        for (j = p_jend; j <= jend_index; j ++)
          {	
	    send_buf[(j - p_jend) * p_numxz + P_2DINDEX(i , k)] = p_vec_val[P_LINDEX_N_GHOST(i, j, k, jbegin_index)];
          }
    
    MPI_Isend(send_buf, send_num, MPI_DOUBLE, mpi_rank + 1, 1, MPI_COMM_WORLD, &request);
  }

  if (mpi_rank != 0) 
  {
    recv_jnum = p_jbegin - jbegin_index + 1;

    recv_num = recv_jnum * p_numxz;

    recv_buf = (double *) malloc(sizeof(double) * recv_num);

    MPI_Recv(recv_buf, recv_num, MPI_DOUBLE, mpi_rank - 1, 1, MPI_COMM_WORLD, &status);

    for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
	for (j = p_jbegin; j <= p_jbegin + ghost_width; j ++)
        {
	  p_vec_val[P_LINDEX_N_GHOST(i,j,k, jbegin_index)] += recv_buf[(j - p_jbegin) * p_numxz + P_2DINDEX(i,k)]; 
        }
  }
  
  if (mpi_rank != mpi_size - 1)
    MPI_Wait(&request, &status);

  if (mpi_rank != 0) free(recv_buf);
  if (mpi_rank != mpi_size - 1) free(send_buf);

   if (mpi_rank != 0) 
  {
    send_jnum = p_jbegin - jbegin_index + 1;

    send_num = send_jnum * p_numxz;

    send_buf = (double *) malloc(sizeof(double) * send_num);
  
    for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
        for (j = jbegin_index; j <= p_jbegin; j ++)
          send_buf[(j - jbegin_index) * p_numxz + P_2DINDEX(i , k)] = p_vec_val[P_LINDEX_N_GHOST(i, j, k, jbegin_index)];
    
    MPI_Isend(send_buf, send_num, MPI_DOUBLE, mpi_rank - 1, 1, MPI_COMM_WORLD, &request);
  }

  if (mpi_rank != mpi_size - 1) 
  {
    recv_jnum = jend_index - p_jend + 1;

    recv_num = recv_jnum * p_numxz;

    recv_buf = (double *) malloc(sizeof(double) * recv_num);

    MPI_Recv(recv_buf, recv_num, MPI_DOUBLE, mpi_rank + 1, 1, MPI_COMM_WORLD, &status);

    int jend_index_begin = p_jend - ghost_width;

    for (i = p_ibegin; i <= p_iend; i ++)
      for (k = p_kbegin; k <= p_kend; k ++)
	      for (j = jend_index_begin; j <= p_jend; j ++)
	          p_vec_val[P_LINDEX_N_GHOST(i,j,k, jbegin_index)] += recv_buf[(j - jend_index_begin) * p_numxz + P_2DINDEX(i,k)]; 
  }
  
  if (mpi_rank != mpi_size - 1)
    MPI_Wait(&request, &status);

  if (mpi_rank != 0) free(send_buf);
  if (mpi_rank != mpi_size - 1) free(recv_buf);

}
