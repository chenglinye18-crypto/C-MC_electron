#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TO3D(i,j,k) ((j) * numx * numz + ((i) * numz + (k)))

int numx, numy, numz, totpoints;
double * lx, * ly, *lz;
double *data, *volume;
char fname[20] = "../data/";

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
  double val[10];
  char buf[20];
  int item_num, index, num;

  printf("input file name:");
  scanf("%s", buf);
  strcat(fname, buf);
  printf("%s\n", fname);
  printf("input data index: ");
  scanf("%d", &index);
  printf("input number of data items per entry:");
  scanf("%d", &item_num);

  fop = fopen(fname, "r");
  for (i = 0;i < totpoints;i ++) {
    fscanf(fop, "%d %d %d", &idx, &idy, &idz);
    for (num = 0; num < item_num ; num ++)
      fscanf(fop, "%lf", &val[num]);
    data[TO3D(idx,idy,idz)] = val[index];
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
  int idx, idy[2], idz[2];
  double range[2];
  FILE * fop;
  char buf[20];

  printf("input x tick value (nm) : ");
  scanf("%lf", &range[0]);
  idx = find_index(lx, numx, range[0]);
  printf("input y range (nm) : ");
  scanf("%lf %lf", &range[0], &range[1]);
  idy[0] = find_index(ly, numy, range[0]);
  idy[1] = find_index(ly, numy, range[1]);
  printf("input z range (nm) : ");
  scanf("%lf %lf", &range[0], &range[1]);
  idz[0] = find_index(lz, numz, range[0]);
  idz[1] = find_index(lz, numz, range[1]);

  printf("input filename for output\n");
  while (getchar() != '\n')
    ;
  scanf("%s", buf);
  printf("output file: %s\n", buf);
  fop = fopen(buf, "w");

  for (i = idy[0]; i <= idy[1]; i ++){
    for (k = idz[0]; k <= idz[1]; k ++)
      fprintf(fop, "%lf %lf %lf\n", ly[i], lz[k], data[TO3D(idx,i,k)]);
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
