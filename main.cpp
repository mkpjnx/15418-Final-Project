#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <getopt.h>
#include "sim.h"
#include "cycleTimer.h"
#include "instrument.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

void write_ppm(grid_t *g, int iter){
  char buf[50];
  sprintf(buf, "out/out%d.ppm", iter);
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

int main(int argc, char** argv){
  bool instrument = false;
  bool verbose = false;
  int steps = 500;
  int runs = 5;
  int dim = 256;
  char opt;
  size_t num_threads = 8;
  while ((opt = getopt(argc, argv, "hvs:r:It:n:")) != -1) {
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
    case 't':
      num_threads = atoi(optarg);
      break;
    case 'n':
      dim = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s [-s steps] \n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  #if defined(_OPENMP)
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(num_threads);
  #endif


  //args: Time steps

  track_activity(instrument); 
  start_activity(ACTIVITY_STARTUP);
  grid_t *g = new_grid(dim,dim);
  initialize_grid(g);
  finish_activity(ACTIVITY_STARTUP);
  double average = 0;
  double start;
  for(int i = 0; i < runs; i ++){
    if (verbose) printf("Run:\t%d\n", i);
    if (instrument) start = CycleTimer::currentSeconds();
    //write_raw(g,i);
    run_grid(g, steps);
    //write_ppm(g, i);
    if (instrument) average += CycleTimer::currentSeconds() - start;
  }
  if (instrument){
    std::cout << "Average time per run:\t" << average/runs << std::endl;
    std::cout << "Average time per step:\t" << average/runs/steps << std::endl;
  }
  show_activity(instrument);
  return 0;
}

