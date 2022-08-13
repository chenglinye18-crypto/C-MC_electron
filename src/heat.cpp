#include "mcmodel.h"

void MeshQuantities::stat_heat_single_step()
{
  int i,j,k;
  double * e_heat, * h_heat, * p_par_charge_val;
  double * e_heat_weight, *h_heat_weight;
  double * vol;

  reduce_add(p_electron_heat, 1);
  reduce_add(p_hole_heat, 1);

/*
  reduce_add(p_e_heat_weight, 1);
  reduce_add(p_h_heat_weight, 1);


  p_electron_heat->ExtractView(&e_heat);
  p_hole_heat->ExtractView(&h_heat);
  p_e_heat_weight->ExtractView(&e_heat_weight);
  p_h_heat_weight->ExtractView(&h_heat_weight);


  for (i = p_ibegin; i <= p_iend; i ++)
    for (k = p_kbegin; k <= p_kend; k ++)
      for (j = p_jbegin_nonoverlap; j <= p_jend_nonoverlap; j ++)
      {
        if (fabs(e_heat_weight[P_LINDEX_ONE_GHOST(i,j,k)]) > 0)
	  e_heat[P_LINDEX_ONE_GHOST(i,j,k)] /= e_heat_weight[P_LINDEX_ONE_GHOST(i,j,k)];
        else {
          if (fabs(e_heat[P_LINDEX_ONE_GHOST(i,j,k)] ) > 0)
            cout << "Warning: electron scattering happened with zero weight\n" 
		 << e_heat[P_LINDEX_ONE_GHOST(i,j,k)] << ' ' << e_heat_weight[P_LINDEX_ONE_GHOST(i,j,k)] << endl;
	}

        if (fabs(h_heat_weight[P_LINDEX_ONE_GHOST(i,j,k)]) > 0)
	  h_heat[P_LINDEX_ONE_GHOST(i,j,k)] /= h_heat_weight[P_LINDEX_ONE_GHOST(i,j,k)];
        else {
	  if (fabs(h_heat[P_LINDEX_ONE_GHOST(i,j,k)] ) > 0)
            cout << "Warning: hole scattering happened with zero weight\n" 
		 << h_heat[P_LINDEX_ONE_GHOST(i,j,k)] << ' ' << h_heat_weight[P_LINDEX_ONE_GHOST(i,j,k)] << endl;

	}
      }
*/
}

void MeshQuantities::heat_density() {

  double *fermi_val, *charge_fac, *qc_fermi_val;
  double *p_par_charge_value;
  double *pot;

  p_charge_fac->ExtractView(&charge_fac);

  p_par_charge->PutScalar(0);
  p_par_charge->ExtractView(&p_par_charge_value);
 
  fill_ghost_cell();

 // particle_to_density(p_par_charge, 2);

  stat_heat_single_step();

  statistic();

  p_electron_heat->PutScalar(0);
  p_hole_heat->PutScalar(0);
/*
  p_e_heat_weight->PutScalar(0);
  p_h_heat_weight->PutScalar(0);
*/

  compute_cell_charge();

  clear_ghost_par();

  if (Flag_SurfaceScatter)
    GetSurfRoughnessPhononScRate();

  t_time->ResetStartTime();
  
  update_particle();

  if (rank == 0)
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

void MeshQuantities::heat_field() {

  pot2field_for_hole(stat_pot_saved);

  pot2field(0, stat_qc_pot_saved);

}

void MeshQuantities::compute_heat() {

  if (rank == 0) {
    cout << endl << "**********************************************" << endl;
    cout << "*" << endl;
    cout << "*" << endl;
    cout << "begin computing heating for fixed potential" << endl;
    cout << "*" << endl;
    cout << "*" << endl;
    cout << "***********************************************" << endl;
  }

  heat_field();

  step = 0;

  flag_heat = true;

  p_electron_heat->PutScalar(0);
  p_hole_heat->PutScalar(0);

  /* simulation loop */
  for (; step < heat_steps; step ++) {

    if (rank == 0) {
      cout << " ===== Heat Computing Process ===== " << endl;
      cout << " ----- Total Heat Step: " << heat_steps << endl;
      cout << "************ step = " << step << " ********" << endl;
      cout << par_num << " particles"  << endl;
    }

    /* self-consistent iteration for one step */
    heat_density();

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

      if (rank == 0)
        cout << "MultipleRefresh Gen : " << tot_mr_gen_num << endl;
    }

    /*particle should be conserved, check it */
   check_par_number();

   /* output statistical result, every stat_step, this may be wrong
    * for the first statistical result, if run from restart */
    if (step % stat_heat_step == (stat_heat_step - 1)){
      output_stat();
      /* output restart result*/
      output_for_restart();
      /*set to zero for the following statistics*/
      set_stat_zero();
    }
  }
}

void MeshQuantities::heat_to_point(double e, Epetra_Vector * p_heat_vec) 
{
    double * p_heat;
    int i = icell, j = jcell, k = kcell;
    double ratio_x, ratio_y, ratio_z;

    ratio_x = (x - lx[i]) / dx[i];
    ratio_y = (y - ly[j]) / dy[j];
    ratio_z = (z - lz[k]) / dz[k];

    if (!(BETWEEN01(ratio_x) && BETWEEN01(ratio_y) && BETWEEN01(ratio_z))) {
      err_message(WRONG_CELL, "particle to density");
//      dump_par_info(*iter);
    }

    p_heat_vec->ExtractView(&p_heat);

    /* distribute the contribution to cell's each node */
    p_heat[P_LINDEX_ONE_GHOST(i,j,k)] += e * (1 - ratio_y) * (1 - ratio_x) * (1 - ratio_z);
    p_heat[P_LINDEX_ONE_GHOST(i + 1, j,k)] += e * (1 - ratio_y) * ratio_x * (1 - ratio_z);
    p_heat[P_LINDEX_ONE_GHOST(i,j,k+1)] += e * (1 - ratio_y) * (1 - ratio_x) * ratio_z;
    p_heat[P_LINDEX_ONE_GHOST(i + 1, j,k + 1)] += e * (1 - ratio_y) * ratio_x * ratio_z;

    p_heat[P_LINDEX_ONE_GHOST(i,j + 1,k)] += e * ratio_y * (1 - ratio_x)* (1 - ratio_z);
    p_heat[P_LINDEX_ONE_GHOST(i + 1, j + 1,k)] += e * ratio_y * ratio_x* (1 - ratio_z);
    p_heat[P_LINDEX_ONE_GHOST(i,j + 1,k+1)] += e * ratio_y * (1 - ratio_x) * ratio_z;
    p_heat[P_LINDEX_ONE_GHOST(i + 1, j + 1,k+1)] += e * ratio_y * ratio_x * ratio_z;
}
