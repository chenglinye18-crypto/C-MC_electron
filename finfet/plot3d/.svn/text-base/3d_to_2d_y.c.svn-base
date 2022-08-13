#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TO3D(i,j,k) ((j) * numx * numz + ((i) * numz + (k)))

int numx, numy, numz, totpoints;
double * lx, * ly, *lz;
double *data, *volume;
char fname[20] = "../data/pot";

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
}

void read_data() {
  int i,j;
  FILE * fop;
  int step, data_idx;
  int idx, idy,idz;
  double val;
  char buf[10];

  printf("input step num: ");
  scanf("%d", &step);
  sprintf(buf, "%d", step); 
  strcat(fname, buf);
  printf("%s\n", fname);

  fop = fopen(fname, "r");
  for (i = 0;i < totpoints;i ++) {
    fscanf(fop, "%d %d %d %lf", &idx, &idy, &idz, &val);
    data[TO3D(idx,idy,idz)] = val;
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

void output() {
  int i,j,k;
  int idx[2], idy, idz[2];
  double range[2];
  FILE * fop;

  printf("input y tick value (nm) : ");
  scanf("%lf", &range[0]);
  idy = find_index(ly, numy, range[0]);
  printf("input x range (nm) : ");
  scanf("%lf %lf", &range[0], &range[1]);
  idx[0] = find_index(lx, numx, range[0]);
  idx[1] = find_index(lx, numx, range[1]);
  printf("input z range (nm) : ");
  scanf("%lf %lf", &range[0], &range[1]);
  idz[0] = find_index(lz, numz, range[0]);
  idz[1] = find_index(lz, numz, range[1]);

  fop = fopen("pot_slice.txt", "w");

  for (i = idx[0]; i <= idx[1]; i ++){
    for (k = idz[0]; k <= idz[1]; k ++)
      fprintf(fop, "%lf %lf %lf\n", lx[i], lz[k], data[TO3D(i, idy, k)]);
    fprintf(fop, "\n");
  }
    
  fclose(fop);
}

int main() {
  int i;

  read_grid();

  read_data();

  output();

  return 0;
}
