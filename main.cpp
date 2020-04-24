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
  int gridsize = 5000;
  char opt;

  int process_count = 1;
  int this_zone = 0;
  int nzone = 0;

  #if MPI
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);
    char* opstring = "hvs:r:Ig:"; //May be unnecessary
  #else
    char* opstring = "hvs:r:Ig:";
  #endif

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
    default:
      fprintf(stderr, "Usage: %s [-s steps] \n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }


  //args: Time steps

  track_activity(instrument); 
  start_activity(ACTIVITY_STARTUP);
  grid_t *g = new_grid(gridsize, gridsize);
  initialize_grid(g);
  state_t *s;
  #if MPI
  if (mpi_master){
      s = divide_grid(g, process_count); //TODO
  } else {
      s = get_divide(process_count, this_zone); //TODO
  }
  init_zone(s);
  #else
    s = (state_t *)malloc(sizeof(state_t *));
  #endif
  s->g = g;
  finish_activity(ACTIVITY_STARTUP);
  double average = 0;
  double start;
  for(int i = 0; i < runs; i ++){
    if (verbose && mpi_master) printf("Run:\t%d\n", i);
    if (instrument && mpi_master) start = CycleTimer::currentSeconds();
    run_grid(s, steps);
    if(mpi_master) write_ppm(g, i);
    if (instrument && mpi_master) average += CycleTimer::currentSeconds() - start;
  }
  if (instrument){
    std::cout << "Average time per run:\t" << average/runs << std::endl;
    std::cout << "Average time per step:\t" << average/runs/steps << std::endl;
  }
  show_activity(instrument);
  #if MPI
    MPI_Finalize();
  #endif
  return 0;
}

