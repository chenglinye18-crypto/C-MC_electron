/*----------------------------------------------------------------------------
every a few dtMultiple Refresh
 details in PhysicsDoc.pdf
----------------------------------------------------------------------------*/


#include "mcmodel.h"

bool MeshQuantities::gen_new_par_mr(vector<Particle> & vec, int new_num, vector<Particle> & new_vec){
  vector<double> weight; 
  int old_num, i, idx;
  double tot_charge, s, aver_charge;
  vector<double>::iterator it;
  double tmp;

  old_num = vec.size();

  if (old_num == 0) return false;

  weight.push_back(0);

  s = 0; 
  tot_charge = 0;
  tmp = 0;
  for (i = 0;i < old_num; i ++){
    s += fabs(vec[i].charge); 
    weight.push_back(s);
    tot_charge += vec[i].charge;
    tmp += vec[i].charge * vec[i].charge;
  }

  if (((old_num > new_num / Rpar) && (old_num < new_num * Rpar)) && (tmp < (tot_charge * tot_charge * Rsqr / old_num)))
    return false;

  for (i = 0;i < new_num; i ++) {
    s = Random() * fabs(tot_charge); 
    it = lower_bound(weight.begin(), weight.end(), s);
    idx = it - weight.begin();

    if ((idx <= old_num) && (idx >= 1)) {
      new_vec.push_back(vec[idx - 1]);
    } else 
      cout << "warning : idx = " << idx << " s = " << s << " old_num = " << old_num << endl
           <<  weight[0] << ' ' << weight[weight.size() - 1] << endl;
  }

  aver_charge = tot_charge / new_num;

  if (new_num != new_vec.size()) {
    cout << "new_num != new_vec.size()" << endl;
    exit(0);
  }

  for (i = 0;i < new_vec.size(); i ++)
    new_vec[i].charge = aver_charge;

  return true;
}

void MeshQuantities::MultipleRefresh(){
  int i, j, k, l;
  int new_ele_num, new_hole_num, before_num;
  int * e_num, * h_num;
  int kill_par;

  list<Particle> * c_par_list;
  list<Particle>::iterator iter;

  vector<Particle> e_vec, h_vec, e_new_vec, h_new_vec;

  c_desired_electron_number->ExtractView(&e_num);
  c_desired_hole_number->ExtractView(&h_num);


  for (i = c_ibegin; i <= c_iend; i ++)
    for (k = c_kbegin; k <= c_kend; k ++)
      for (j = c_jbegin; j <= c_jend; j ++)
        {
	  /* get the cell's particle list */	
	  c_par_list = &par_list[C_LINDEX_GHOST_ONE(i,j,k)];

	  e_new_vec.resize(0);
	  h_new_vec.resize(0);
	  e_vec.resize(0);
	  h_vec.resize(0);
	  kill_par = 0;

	  /*for each particle */
	  for (iter = c_par_list->begin(); iter != c_par_list->end(); iter ++)
            if (fabs(iter->charge) >= pckill){
	      if (iter->par_type == PELEC) 
		e_vec.push_back(*iter);
              else 
		h_vec.push_back(*iter);
            } else kill_par ++;


          new_ele_num = e_num[C_LINDEX_GHOST_ONE(i,j,k)]; 
          new_hole_num = h_num[C_LINDEX_GHOST_ONE(i,j,k)]; 

	  before_num = c_par_list->size();

	  c_par_list->clear();

	  if (gen_new_par_mr(e_vec, new_ele_num, e_new_vec)) {
	      for (l = 0; l < e_new_vec.size(); l ++)
		c_par_list->push_back(e_new_vec[l]);
	    } else {
	      for (l = 0; l < e_vec.size(); l ++)
		c_par_list->push_back(e_vec[l]);
	    }

          if (gen_new_par_mr(h_vec, new_hole_num, h_new_vec)) {
	      for (l = 0; l < h_new_vec.size(); l ++)
		c_par_list->push_back(h_new_vec[l]);
	    } else {
	      for (l = 0; l < h_vec.size(); l ++)
		c_par_list->push_back(h_vec[l]);
	    }

	  mr_gen_num += (c_par_list->size() - before_num);
	}
}
