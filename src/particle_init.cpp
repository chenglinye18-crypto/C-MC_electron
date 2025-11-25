/**
 * @file particle_init.cpp
 * @brief 粒子初始化相关实现，从 mcmodel.cpp 拆分便于单独阅读
 */

/**
 * @brief 计算每个 cell 初始需要的电子/空穴粒子数量
 *
 * 输入来源：
 *  - c_init_electron_charge / c_init_hole_charge：根据掺杂得到的电荷分布
 *  - sum_charge[PELEC/PHOLE] 与 electron_number/hole_number：目标全局粒子数
 *  - c_desired_* 与 default_*：用户指定或默认的最低粒子数
 *
 * 输出结果：
 *  - 写回 c_electron_num / c_hole_num，后续 init_particle_data 会按此生成粒子
 */
void MeshQuantities::compute_par_num() {

   int * material;
   int i,j,k;
   Particle newpar;
   double lamda, par_charge;
   int pnum;
   int iptype, ipar;
   int * cell_electron_num, * cell_hole_num, * desired_ele_num, * desired_hole_num;
   double * electron_charge;
   double * hole_charge;

   c_material->ExtractView(&material);    //每个cell的材料类型
   c_electron_num->ExtractView(&cell_electron_num);  //输出每个cell中电子数量
   c_hole_num->ExtractView(&cell_hole_num);    //输出每个cell中空穴数量

   // 这两者的初始化在 init_cell_data 中完成
   c_desired_electron_number->ExtractView(&desired_ele_num);  //最小粒子数要求
   c_desired_hole_number->ExtractView(&desired_hole_num);  //最小空穴要求

   /*initialized in init_cell_data(), according to the dopping density of
    * donors and acceptors.
    * */
   c_init_electron_charge->ExtractView(&electron_charge);
   c_init_hole_charge->ExtractView(&hole_charge);

   /*compute particle number for each cell */
   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++)
       /* only silicon cell contains carriers */
	 if(material[C_LINDEX_GHOST_ONE(i,j,k)] == SILICON) {
       /* electron_number and hole_number are set by the config files as the
	* total number of electrons and holes to be simulated*/
         /* 按电荷占比分配粒子数，再 +1 保证至少一个 */
         cell_electron_num[C_LINDEX_GHOST_ONE(i,j,k)] =
           int (1e-10 + electron_charge[C_LINDEX_GHOST_ONE(i,j,k)] / sum_charge[PELEC] * electron_number) + 1; 
         cell_hole_num[C_LINDEX_GHOST_ONE(i,j,k)] =
           int (1e-10 + hole_charge[C_LINDEX_GHOST_ONE(i,j,k)] / sum_charge[PHOLE] * hole_number) + 1;

         /* 如果自动计算的数量超过了用户下限，则提升下限，避免后面被覆盖掉 */
         if (cell_electron_num[C_LINDEX_GHOST_ONE(i,j,k)] > desired_ele_num[C_LINDEX_GHOST_ONE(i,j,k)])
	   desired_ele_num[C_LINDEX_GHOST_ONE(i,j,k)] = cell_electron_num[C_LINDEX_GHOST_ONE(i,j,k)];

         if (cell_hole_num[C_LINDEX_GHOST_ONE(i,j,k)] > desired_hole_num[C_LINDEX_GHOST_ONE(i,j,k)])
	   desired_hole_num[C_LINDEX_GHOST_ONE(i,j,k)] = cell_hole_num[C_LINDEX_GHOST_ONE(i,j,k)];

	 /* 至少保证默认的粒子数，下限保护 */
	 if (desired_ele_num[C_LINDEX_GHOST_ONE(i,j,k)] < default_electron_num)
	   desired_ele_num[C_LINDEX_GHOST_ONE(i,j,k)] = default_electron_num;

	 if (desired_hole_num[C_LINDEX_GHOST_ONE(i,j,k)] < default_hole_num)
	   desired_hole_num[C_LINDEX_GHOST_ONE(i,j,k)] = default_hole_num;

       /* count the total carrier number on our local process, excluding the
	* ghost cells */
       }
 }

/**
 * @brief 统计本进程及全局的粒子总数，用于日志/调试
 */
void MeshQuantities::compute_total_par_num(){
   int i,j,k;
   int * cell_electron_num, * cell_hole_num;

   local_par_num = 0;

   c_electron_num->ExtractView(&cell_electron_num);
   c_hole_num->ExtractView(&cell_hole_num);

   for (i = c_ibegin; i <= c_iend; i ++)
     for (k = c_kbegin; k <= c_kend; k ++)
       for (j = c_jbegin; j <= c_jend; j ++)
	 local_par_num += cell_electron_num[C_LINDEX_GHOST_ONE(i,j,k)] + cell_hole_num[C_LINDEX_GHOST_ONE(i,j,k)];

   /*count the total number by mpi_allreduce */
   MPI_Allreduce(&local_par_num, &par_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if (mpi_rank == 0)
     cout << "Particle Number : " << par_num << endl;

 }

/** 
 * @brief initialize particles array
 *
 * 处理流程：
 * 1) 通过 compute_par_num() 计算各 cell 需要的粒子数量；
 * 2) 若是重启，读取保存的电荷分布，以维持原有的局部中性；
 * 3) 对每个硅 cell，根据需要的粒子数在单元内随机撒点，设置电荷、动量状态并放入 par_list。
 *
 * 生成结果：
 * par_list[ cell_idx ] 中存放该 cell 的所有初始粒子，每个粒子带有位置、索引、charge、k-state。
 */
 void MeshQuantities::init_particle_data() {
   int * material;
   int i,j,k;
   Particle newpar;
   double lamda, par_charge;
   int pnum;
   int iptype, ipar;
   int * cell_electron_num, * cell_hole_num; 
   double * electron_charge;
   double * hole_charge;
   double * ele_charge_value, * hole_charge_value;
   int * cell_contact;

   Epetra_Vector c_ele_charge_init(*c_map_ghost);
   Epetra_Vector c_hole_charge_init(*c_map_ghost);

   compute_par_num();  // 计算每个cell中的电子和空穴数量

   if (Flag_restart) {
     read_par_info_for_restart(&c_ele_charge_init, c_electron_num,&c_hole_charge_init, c_hole_num);
     c_ele_charge_init.ExtractView(&ele_charge_value);
     c_hole_charge_init.ExtractView(&hole_charge_value);
   }

   //计算每个cell中的电子和空穴数量总和
   compute_total_par_num();

   c_material->ExtractView(&material);

   c_electron_num->ExtractView(&cell_electron_num);
   c_hole_num->ExtractView(&cell_hole_num);

   c_init_electron_charge->ExtractView(&electron_charge);
   c_init_hole_charge->ExtractView(&hole_charge);

   c_attached_contact->ExtractView(&cell_contact);

     /* for each cells and each particle type */
   for(iptype = 0; iptype < 2; iptype ++)
     for (i = c_ibegin; i <= c_iend; i ++)
       for (k = c_kbegin; k <= c_kend; k ++)
	 for (j = c_jbegin; j <= c_jend; j ++)
	 /* silicon cells only*/
         if (material[C_LINDEX_GHOST_ONE(i,j,k)] == SILICON) 
         {
	   /* 非重启或连接接触的 cell：用当前中性电荷分配给每个粒子 */
	   if ((!Flag_restart) || (cell_contact[C_LINDEX_GHOST_ONE(i,j,k)] > 0)) {
	     if (iptype == 0){/*particles share the same charge*/
	       pnum = cell_electron_num[C_LINDEX_GHOST_ONE(i,j,k)];
	       par_charge = electron_charge[C_LINDEX_GHOST_ONE(i,j,k)] / pnum;
	     }
	     else{/*holes share the same charge */
	       pnum = cell_hole_num[C_LINDEX_GHOST_ONE(i,j,k)];
	       par_charge = hole_charge[C_LINDEX_GHOST_ONE(i,j,k)] / pnum;
	     }
	   } else {
	     /* 重启且非接触区：沿用旧电荷分布 */
	     if (iptype == 0){/*particles share the same charge*/
	       pnum = cell_electron_num[C_LINDEX_GHOST_ONE(i,j,k)];
	       if (pnum == 0) continue;
	       par_charge = - ele_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] / pnum;
	     } else{/*holes share the same charge */
	       pnum = cell_hole_num[C_LINDEX_GHOST_ONE(i,j,k)];
	       if (pnum == 0) continue;
	       par_charge = hole_charge_value[C_LINDEX_GHOST_ONE(i,j,k)] / pnum;
	     }
	   }

	   for (ipar = 0; ipar < pnum; ipar ++){ 

             newpar.par_id = ipar;

             lamda = Random();
	     if (!BETWEEN01(lamda)){
	       err_message(WRONG_RANDOM, "in init_particle");
	       exit(0);
	     }

             /* 在当前 cell 内沿 x/y/z 随机均匀撒点 */
             newpar.x = lamda * lx[i] + (1 - lamda) * lx[i + 1];
             lamda = Random();
	     if (!BETWEEN01(lamda)){
	       err_message(WRONG_RANDOM, "in init_particle");
	       exit(0);
	     }
             newpar.y = lamda * ly[j] + (1 - lamda) * ly[j + 1];
             lamda = Random();

	     if (!BETWEEN01(lamda)){
	       err_message(WRONG_RANDOM, "in init_particle");
	       exit(0);
	     }

             newpar.z = lamda * lz[k] + (1 - lamda) * lz[k + 1];

             newpar.charge = par_charge;

             newpar.i = i;
             newpar.j = j;
             newpar.k = k;
             newpar.par_type = iptype;
             newpar.seed = seed;
             newpar.left_time = dt;
	     select_kstate(&newpar, 0);
             par_list[C_LINDEX_GHOST_ONE(i,j,k)].push_back(newpar);
	 }
       }
 }
