#include "utils.h"

int numx, numy, numz, totpoints;
double * lx, * ly, *lz;
double *data, *volume;

char conf_name[MAXLINE + 1] = "config.txt";
char datafile_name[MAXLINE + 1];
char output_name[MAXLINE + 1];
char gridfile_name[MAXLINE + 1];
char fig_name[MAXLINE + 1];

int entry_num, entry_idx;
int section_type;

double begin_range[2], end_range[2], pos;
int begin_index[2], end_index[2], pos_index;

void read_config(FILE * conf_fd) {
  fscanf(conf_fd, "%s", datafile_name);
  fscanf(conf_fd, "%d %d", &entry_num, &entry_idx);
  fscanf(conf_fd, "%d", &section_type);
  fscanf(conf_fd, "%lf %lf %lf %lf %lf", &pos, &begin_range[0], &end_range[0], &begin_range[1], &end_range[1]);
  fscanf(conf_fd, "%s", output_name);
  fscanf(conf_fd, "%s", fig_name);

#ifdef DEBUG
  printf("datafile_name : %s\n entry_num = %d\n entry_idx = %d\n section_type = %d \n output_name = %s\n", \
      datafile_name, entry_num, entry_idx, section_type, output_name);
#endif

}

void read_grid() {
  FILE * fop;
  int i;

  fop = fopen(gridfile_name, "r");
  if (fop == NULL)
    err_sys("can not open grid file");

  fscanf(fop, "%d", &numx);
  lx = (double *) malloc(numx * sizeof(double));

  for (i = 0;i < numx; i ++){
    fscanf(fop, "%lf", &lx[i]);
    lx[i] *= 1000;
  }

  fscanf(fop, "%d", &numy);
  ly = (double *) malloc(numy * sizeof(double));

  for (i = 0;i < numy; i ++){
    fscanf(fop, "%lf", &ly[i]);
    ly[i] *= 1000;
  }

  fscanf(fop, "%d", &numz);
  printf("%d %d %d\n", numx, numy, numz);
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
  int idx, idy, idz;
  double val[MAX_ENTRY_NUM];
  int num;

  fop = fopen(datafile_name, "r");
  if (fop == NULL)
    err_sys("can not open datafile file %s", datafile_name);

  for (i = 0;i < totpoints; i ++) {
    fscanf(fop, "%d %d %d", &idx, &idy, &idz);
    for (num = 0; num < entry_num; num ++)
      fscanf(fop, "%lf", &val[num]);
    data[TO3D(idx,idy,idz)] = val[entry_idx];
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

  fop = fopen(output_name, "w");

  switch (section_type) {
case DIR_X:
    pos_index = find_index(lx, numx, pos);
    begin_index[0] = find_index(ly, numy, begin_range[0]);
    end_index[0] = find_index(ly, numy, end_range[0]);
    begin_index[1] = find_index(lz, numz, begin_range[1]);
    end_index[1] = find_index(lz, numz, end_range[1]);

#ifdef DEBUG
    printf("x pos_index = %d, y dir from %d to %d , z dir from %d to %d\n", pos_index, begin_index[0], end_index[0], begin_index[1], end_index[1]);
#endif

    for (j = begin_index[0]; j <= end_index[0]; j ++){
      for (k = begin_index[1]; k <= end_index[1]; k ++)
	fprintf(fop, "%lf %lf %lf\n", ly[j], lz[k], data[TO3D(pos_index, j, k)]);
      fprintf(fop, "\n");
    }
    break;

case DIR_Y: 
    pos_index = find_index(ly, numy, pos);
    begin_index[0] = find_index(lx, numx, begin_range[0]);
    end_index[0] = find_index(lx, numx, end_range[0]);
    begin_index[1] = find_index(lz, numz, begin_range[1]);
    end_index[1] = find_index(lz, numz, end_range[1]);

#ifdef DEBUG
    printf("y pos_index = %d, x dir from %d to %d , z dir from %d to %d\n", pos_index, begin_index[0], end_index[0], begin_index[1], end_index[1]);
#endif

    for (i = begin_index[0]; i <= end_index[0]; i ++){
      for (k = begin_index[1]; k <= end_index[1]; k ++)
	fprintf(fop, "%lf %lf %lf\n", lx[i], lz[k], data[TO3D(i, pos_index, k)]);
      fprintf(fop, "\n");
    }
    break;
case DIR_Z: 
    pos_index = find_index(lz, numz, pos);
    begin_index[0] = find_index(lx, numx, begin_range[0]);
    end_index[0] = find_index(lx, numx, end_range[0]);
    begin_index[1] = find_index(ly, numy, begin_range[1]);
    end_index[1] = find_index(ly, numy, end_range[1]);

#ifdef DEBUG
    printf("z pos_index = %d, x dir from %d to %d , y dir from %d to %d\n", pos_index, begin_index[0], end_index[0], begin_index[1], end_index[1]);
#endif

    for (i = begin_index[0]; i <= end_index[0]; i ++){
      for (j = begin_index[1]; j <= end_index[1]; j ++)
	fprintf(fop, "%lf %lf %lf\n", lx[i], ly[j], data[TO3D(i, j, pos_index)]);
      fprintf(fop, "\n");
    }
    break;
default: 
      err_sys("wrong section type");
      break;
  }
   
  fclose(fop);
}

int main(int argc, char ** argv) {
  int i;
  FILE * conf_fd;
  int num_fig;

  if (argc == 2)
    strncpy(conf_name, argv[1], sizeof(argv[1])); 

  conf_fd = fopen(conf_name, "r");

  if (conf_fd == NULL)
    err_sys("can not open config file %s", conf_name);

  fscanf(conf_fd, "%s", gridfile_name);
  fscanf(conf_fd, "%d", &num_fig);

  read_grid();

  for (i = 0;i < num_fig; i ++) {

    read_config(conf_fd);

    read_data();

    output();

  }

  return 0;
}
