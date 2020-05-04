#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <getopt.h>
#include "sim.h"
#include "setup.h"
#include "cycleTimer.h"
#include "instrument.h"

void Exit(){
  #if MPI
    MPI_Finalize();
  #endif
  exit(0);
}

void usage(){
  printf("-h\t help\n");
  printf("-v\t verbose mode\n");
  printf("-s S\t steps per run\n");
  printf("-r R\t runs\n");
  printf("-I\t verbose mode\n");
  printf("-g G\t grid size\n");
  printf("-d D\t horizantal divisions\n");
}

void write_raw(grid_t *g, int iter){
  char buf[50];
  sprintf(buf, "out/out%d.txt", iter);
  FILE *fp = fopen(buf, "wb");

  if (!fp) {
      fprintf(stderr, "Error: could not open file for write\n");
      exit(1);
  }
  for (int j=g->nrow-1; j>=0; j--) {
    for (int i=0; i< g->ncol; i++) {

      double dub = g->v[GINDEX(g,j,i)];
      fprintf(fp, "%d,%d: %a\n", j, i, dub);
    }
  }
    fclose(fp);
}

void write_ppm(grid_t *g, int iter){
  char buf[50];
  #if MPI
    sprintf(buf, "out/MPI%d.ppm", iter);
  #else
    sprintf(buf, "out/out%d.ppm", iter);
  #endif
  FILE *fp = fopen(buf, "wb");

  if (!fp) {
      fprintf(stderr, "Error: could not open file for write\n");
      exit(1);
  }
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", g->ncol, g->nrow);
  fprintf(fp, "255\n");
  for (int j=g->nrow-1; j>=0; j--) {
    for (int i=0; i< g->ncol; i++) {

      double dub = g->v[GINDEX(g,j,i)] * 3;

      char val[3];
      val[0] = static_cast<char>(255. * CLAMP(dub, 0., 1.));
      val[1] = static_cast<char>(255. * CLAMP(dub, 0., 1.));
      val[2] = static_cast<char>(255. * CLAMP(dub, 0., 1.));

      fputc(val[0], fp);
      fputc(val[1], fp);
      fputc(val[2], fp);
    }
  }
    fclose(fp);
}

int main(int argc, char** argv){
  bool instrument = false;
  bool verbose = false;
  int steps = 500;
  int runs = 5;
  int gridsize = 500;
  char opt;

  int process_count = 1;
  int horizantal_divisions = 1;
  int this_zone = 0;
  int nzone = 0;

  #if MPI
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);
  #endif
  char* opstring = "hvs:r:Ig:d:";

  bool mpi_master = this_zone == 0;

  
  while ((opt = getopt(argc, argv, opstring)) != -1) {
    switch (opt) {
    case 'h':
      printf("h for howdy\n");
      break;
    case 'v':
      verbose = true;
      break;
    case 's':
      steps = atoi(optarg);
      break;
    case 'r':
      runs = atoi(optarg);
      break;
    case 'I':
      instrument = true;
      break;
    case 'g':
      gridsize = atoi(optarg);
      break;
    case 'd':
      horizantal_divisions = atoi(optarg);
      break;
    default:
      printf("Unknown argument");
      usage();
      Exit();
    }
  }


  //args: Time steps

  track_activity(instrument); 
  if (mpi_master) start_activity(ACTIVITY_STARTUP);
  state_t *s = init_zone(gridsize, gridsize, process_count, this_zone, horizantal_divisions);
  if (s == NULL) {
    fprintf(stderr, "Zone not initialized\n");
    Exit();
  }
  initialize_grid(s);
  if (mpi_master) finish_activity(ACTIVITY_STARTUP);
  double average = 0;
  double start;
  grid_t *new_g;
  for(int i = 0; i < runs; i ++){
    if (verbose && mpi_master) printf("Run:\t%d\n", i);
    new_g = run_grid(s, steps, M_REDBLACK);
    if(mpi_master) write_ppm(new_g, i);
  }
  if (mpi_master){
    show_activity(instrument);
  }
  #if MPI
    MPI_Finalize();
  #endif
  return 0;
}

