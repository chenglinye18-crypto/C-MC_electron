#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TO3D(i,j,k) ((j) * numx * numz + ((i) * numz + (k)))

int numx, numy, numz, totpoints;
double * lx, * ly, *lz;
double *data, *volume, *dens;
char fname[50] = "../data/Electron";
int range_idx[2], range_idy[2], range_idz[2];
double range[2];
char output_file[50];
int step;
int data_idx;

void read_grid() {
  FILE * fop;
  int i;
  fop = fopen("../lgrid.txt", "r");

  fscanf(fop, "%d", &numx);
  lx = (double *) malloc(numx * sizeof(double));

  for (i = 0;i < numx; i ++){
    fscanf(fop, "%lf", &lx[i]);
    lx[i] *= 1000;
  }

  fscanf(fop, "%d", &numy);
  printf("%d\n", numy);
  ly = (double *) malloc(numy * sizeof(double));

  for (i = 0;i < numy; i ++){
    fscanf(fop, "%lf", &ly[i]);
    ly[i] *= 1000;
  }

  fscanf(fop, "%d", &numz);
  printf("%d\n", numz);
  lz = (double *)malloc(numz * sizeof(double));

  for (i = 0;i < numz; i ++){
    fscanf(fop, "%lf", &lz[i]);
    lz[i] *= 1000;
  }

  fclose(fop);

  totpoints = numx * numy * numz;
  volume = (double *) malloc(totpoints * sizeof(double));
  data = (double *) malloc(totpoints * sizeof(double));
  dens = (double *) malloc(totpoints * sizeof(double));
}

void read_volume() {
  int i;
  FILE * fop;
  int idx, idy,idz;
  double val;

  fop = fopen("../data/pvolume", "r");

  for (i = 0;i < totpoints;i ++) {
    fscanf(fop, "%d %d %d %lf", &idx, &idy, &idz, &val);
    volume[TO3D(idx,idy,idz)] = val;
    //if (val > 0) printf("%d %d %d %lf\n", idx, idy, idz, val);
  }
  fclose(fop);
}

void read_data() {
  int i,j;
  FILE * fop;
  int idx, idy,idz;
  double val[8];
  char buf[30];

  fop = fopen("int_input.txt","r");
  fscanf(fop, "%d", &step);
  fscanf(fop, "%d", &data_idx);
  sprintf(buf, "%d", step); 
  strcat(fname, buf);
  printf("%s\n", fname);
  fscanf(fop, "%s", output_file);

  fscanf(fop, "%lf %lf", &range[0], &range[1]);
  range_idy[0] = find_index(ly, numy, range[0]);
  range_idy[1] = find_index(ly, numy, range[1]);

  fscanf(fop, "%lf %lf", &range[0], &range[1]);
  range_idx[0] = find_index(lx, numx, range[0]);
  range_idx[1] = find_index(lx, numx, range[1]);

  fscanf(fop, "%lf %lf", &range[0], &range[1]);
  range_idz[0] = find_index(lz, numz, range[0]);
  range_idz[1] = find_index(lz, numz, range[1]);

  fclose(fop);

  fop = fopen(fname, "r");
  for (i = 0;i < totpoints;i ++) {
    fscanf(fop, "%d %d %d", &idx, &idy, &idz);
    for (j = 0;j < 5; j ++)
      fscanf(fop, "%lf", &val[j]);
    data[TO3D(idx,idy,idz)] = val[data_idx];
    if (data_idx == 1) 
      data[TO3D(idx, idy, idz)] *= 100;
    dens[TO3D(idx,idy,idz)] = fabs(val[4]);
  }
  fclose(fop);
}

int find_index(double * l, int num, double val) {
  int i;
  for (i = 0;i < num; i ++){
    if (l[i] - val >= -1e-9) return i;
  }
  return num - 1;
}

void integral() {
  int i,j,k;
  FILE * fop;
  double * weight, * int_data;

  weight = (double *) malloc(numy * sizeof(double));
  int_data = (double *) malloc(numy * sizeof(double));

  for (j = range_idy[0]; j <= range_idy[1]; j ++) {
    weight[j] = 0;
    int_data[j] = 0;
    for (i = range_idx[0]; i <= range_idx[1]; i ++)
      for (k = range_idz[0]; k <= range_idz[1]; k ++){
	if (data_idx == 4){
	  int_data[j] += data[TO3D(i,j,k)] * volume[TO3D(i,j,k)];
	  weight[j] += volume[TO3D(i,j,k)];
	} else {
	    int_data[j] += data[TO3D(i,j,k)] * volume[TO3D(i,j,k)] * dens[TO3D(i,j,k)];
	    weight[j] += volume[TO3D(i,j,k)] * dens[TO3D(i,j,k)];
	}
      }
    int_data[j] /= weight[j];
  }

  fop = fopen(output_file, "w");
  for (i = range_idy[0]; i <= range_idy[1]; i ++) 
    fprintf(fop, "%lf %lf\n", ly[i], int_data[i]);
  fclose(fop);

}

int main() {
  int i;

  read_grid();

  read_volume();

  read_data();

  integral();

  return 0;
}
