#include "mcmodel.h"

void MeshQuantities::init_vector_entry(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK,int EndK,
                                       double entry_value, Epetra_Vector * c_vector) 
{
  int i,j,k;
  double * vector_array;
  
  c_vector->ExtractView(&vector_array);
    
  if ((BeginI >= 0) && (EndI <= p_numx - 1) && (BeginI < EndI)
   && (BeginK >= 0) && (EndK <= p_numz - 1) && (BeginK < EndK)
   && (BeginJ >= 0) && (EndJ <= p_numy - 1) && (BeginJ < EndJ)) {

    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);
    
    for(i = BeginI; i < EndI; i++)
      for(k = BeginK; k < EndK; k++)
	for(j = BeginJ; j < EndJ; j++)
	  if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
	    vector_array[C_LINDEX_GHOST_ONE(i,j,k)] = entry_value;
  }  else 
      cout << "Set Vector : Overflow or ijk error" << endl;
}


void MeshQuantities::init_vector_entry(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK,int EndK,
                                       int entry_value, Epetra_IntVector * c_vector) 
{
  int i,j,k;
  int * vector_array;

  c_vector->ExtractView(&vector_array);
    
  if ((BeginI >= 0) && (EndI <= p_numx - 1) && (BeginI < EndI)
   && (BeginK >= 0) && (EndK <= p_numz - 1) && (BeginK < EndK)
   && (BeginJ >= 0) && (EndJ <= p_numy - 1) && (BeginJ < EndJ)) {

    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);
    
    for(i = BeginI; i < EndI; i++)
      for(k = BeginK; k < EndK; k++)
	for(j = BeginJ; j < EndJ; j++)
	  if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
	    vector_array[C_LINDEX_GHOST_ONE(i,j,k)] = entry_value;
	  
  }
  else 
    cout << "Set int vector : Overflow or ij error" << endl;
}

void MeshQuantities::SetMRCube(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK,int EndK,
                        int mr0,int mr1,int mr2,int mr3, int mr4, int mr5)
{
  int i,j,k;

  int * mr_value;
  
  c_motion_rule->ExtractView(&mr_value);
                         
  if ((BeginI >= 0) && (EndI <= p_numx - 1) && (BeginI < EndI)
   && (BeginK >= 0) && (EndK <= p_numz - 1) && (BeginK < EndK)
   && (BeginJ >= 0) && (EndJ <= p_numy - 1) && (BeginJ < EndJ)) {

    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);
    
    for(i = BeginI; i < EndI; i++)
      for(k = BeginK; k < EndK; k++)
	for(j = BeginJ; j < EndJ; j++)
	  if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
          {
            mr_value[C_LINDEX_GHOST(i, j, k,0, 6)] = mr0;
            mr_value[C_LINDEX_GHOST(i, j,k, 1, 6)] = mr1;
            mr_value[C_LINDEX_GHOST(i, j,k, 2, 6)] = mr2;
            mr_value[C_LINDEX_GHOST(i, j,k, 3, 6)] = mr3;
            mr_value[C_LINDEX_GHOST(i, j,k, 4, 6)] = mr4;
            mr_value[C_LINDEX_GHOST(i, j,k, 5, 6)] = mr5;
          }
  }
  else 
    cout<<"SetMRCube : Overflow or ijk error"<<endl;
}

void MeshQuantities::SetMRPlane(int BeginI,int EndI,int BeginJ,int EndJ,int BeginK,int EndK,

                         int idn,int mr)
{
  int * mr_value;
  int i,j,k;
  
  c_motion_rule->ExtractView(&mr_value);
                           
  if ((BeginI >= 0) && (EndI <= p_numx - 1) && (BeginI <=EndI)
   && (BeginK >= 0) && (EndK <= p_numz - 1) && (BeginK <=EndK)
   && (BeginJ >= 0) && (EndJ <= p_numy - 1) && (BeginJ <=EndJ)) {

    if (BeginI==EndI) 
    {
      if ((BeginJ==EndJ)|| (BeginK == EndK))
        cout << "SetMRPlane : Must be a Plane";
      else {
        if (idn==UP) {
	  if (BeginI <= c_numx - 1) {
	    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
	    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);

	    for(k = BeginK; k < EndK; k++)
	      for(j = BeginJ; j < EndJ; j++)
		if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
                      mr_value[C_LINDEX_GHOST(BeginI, j, k, UP, 6)] = mr;
	  } else
	    cout << "SetMRPlane: Cannot set UP montion rule of a UP surface" << endl;
	} else if (idn == DOWN) {
	  if (BeginI > 0) {
	    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
	    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);

	    for(k = BeginK; k < EndK; k++)
	      for(j = BeginJ; j < EndJ; j++)
		if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
		  mr_value[C_LINDEX_GHOST(BeginI - 1, j, k, DOWN, 6)] = mr;
	    } else
	    cout << " SetMRPlane: Cannot set DOWM montion rule of a DOWN surface"<<endl;
	} else
	  cout << "SetMRPlane : Wrong Direction";
      }
    } else 
      if(BeginJ == EndJ) {
	if (BeginK == EndK)
	  cout << "SetMRPlane : Must be a Plane";
	else {
	if (idn == LEFT) {
          if ((BeginJ >= c_jbegin_ghost) && (BeginJ <= c_jend_ghost)) {
	    for(k = BeginK; k < EndK; k++)
	      for(i = BeginI;i < EndI; i++)
		mr_value[C_LINDEX_GHOST(i, BeginJ,k, LEFT, 6)] = mr;
	  }
	} else if (idn == RIGHT) {
	  if ((BeginJ - 1>= c_jbegin_ghost) && (BeginJ - 1  <= c_jend_ghost)) {
	    for(k = BeginK; k < EndK; k++)
	      for(i = BeginI; i< EndI; i++)
		mr_value[C_LINDEX_GHOST(i, BeginJ-1, k, RIGHT, 6)] = mr;
	  }
	} else 
            cout<<"SetMRPlane : Wrong Direction";
	}
      }
    else
    if (BeginK==EndK) {
      if (idn== BACK) {
	if (BeginK > 0) {
	    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
	    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);

	    for(i = BeginI; i < EndI; i++)
	      for(j = BeginJ; j < EndJ; j++)
		if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
                      mr_value[C_LINDEX_GHOST(i, j, BeginK - 1, BACK, 6)] = mr;
	  } else
	    cout << "SetMRPlane: Cannot set BACK montion rule of a BACK surface" << endl;
	} else if (idn == FRONT) {
	  if (BeginK < c_numz - 1) {
	    BeginJ = P_LOCAL_LRANGE_GHOST(BeginJ);
	    EndJ = P_LOCAL_RRANGE_GHOST(EndJ);

	    for(i = BeginI; i < EndI; i++)
	      for(j = BeginJ; j < EndJ; j++)
		if ((j >= c_jbegin_ghost) && (j <= c_jend_ghost))
		  mr_value[C_LINDEX_GHOST(i, j, BeginK, FRONT, 6)] = mr;
	    } else
	    cout << " SetMRPlane: Cannot set FRONT montion rule of a FRONT surface"<<endl;
	} else
	  cout << "SetMRPlane : Wrong Direction";
      }

  }
}



