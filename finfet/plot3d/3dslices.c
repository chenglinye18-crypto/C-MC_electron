#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TO3D(i,j,k) ((j) * numx * numz + ((i) * numz + (k)))

int numx, numy, numz, totpoints;
double * lx, * ly, *lz;
double *data, *volume;
char fname[20] = "../data/";
char out_fname[20];
int item_num, data_index, num;
int num_slices;
int idx[2], idz[2];
double range[2];
int idy[30];


void read_input() {
  FILE * fop;
  char buf[20];
  int i;

  //printf("input file name:");
  fop = fopen("slices_input.txt","r");
  fscanf(fop, "%s", buf);
  strcat(fname, buf);
  printf("%s\n", fname);
  //printf("input data index: ");
  fscanf(fop, "%d %d", &data_index, &item_num);
  //printf("input number of data items per entry:");
  
  fscanf(fop, "%d", &num_slices);

  //printf("input y tick values (nm) : ");
  for (i = 0;i < num_slices;i ++) {
    fscanf(fop, "%lf", &range[0]);
    idy[i] = find_index(ly, numy, range[0]);
  }

  //printf("input x range (nm) : ");
  fscanf(fop, "%lf %lf", &range[0], &range[1]);
  idx[0] = find_index(lx, numx, range[0]);
  idx[1] = find_index(lx, numx, range[1]);

  //printf("input z range (nm) : ");
  fscanf(fop, "%lf %lf", &range[0], &range[1]);
  idz[0] = find_index(lz, numz, range[0]);
  idz[1] = find_index(lz, numz, range[1]);

 //printf("input filename for output\n");
 //  while (getchar() != '\n')
 //   ;
  fscanf(fop, "%s", out_fname);
  printf("output file: %s\n", out_fname);
}

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
  int step;
  int idx, idy,idz;
  double val[10];
  char buf[20];

  fop = fopen(fname, "r");
  for (i = 0;i < totpoints;i ++) {
    fscanf(fop, "%d %d %d", &idx, &idy, &idz);
    for (num = 0; num < item_num ; num ++) 
      fscanf(fop, "%lf", &val[num]);
    data[TO3D(idx,idy,idz)] = val[data_index];
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
  FILE * fop;
  char out_slice_fname[30], snum[5];

  //printf("input num of slices: ");

  for (j = 0;j < num_slices; j ++){
    strcpy(out_slice_fname, out_fname);
    sprintf(snum, "%d", j);
    strcat(out_slice_fname, snum);
    strcat(out_slice_fname, ".txt");
    printf("%s\n", out_slice_fname);
    fop = fopen(out_slice_fname, "w");
    for (i = idx[0]; i <= idx[1]; i ++) {
      for (k = idz[0]; k <= idz[1]; k ++)
	fprintf(fop, "%lf %lf %lf %lf\n", lz[k],ly[idy[j]], lx[i], data[TO3D(i,idy[j], k)]);
      fprintf(fop, "\n");
    }
    fclose(fop);
  }
}

int main() {
  int i;

  read_grid();

  read_input();

  read_data();

  output();

  return 0;
}
