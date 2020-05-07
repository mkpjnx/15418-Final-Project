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

void usage(){
  printf("-h\t help\n");
  printf("-v\t verbose mode\n");
  printf("-s S\t steps per run \t\t Default: 500\n");
  printf("-I\t instrument mode \n");
  printf("-g G\t grid size \t\t Default: 256\n");
  printf("-t T\t number of threads \t Default: 8\n");
  printf("-u U\t Update function 0 for Jacobi 1 for Red Black \t Default: 0\n");
}


int main(int argc, char** argv){
  //For performance modeling
  bool instrument = false, verbose = false;

  //Default values for the grid and simulation paramaters
  int steps = 500, runs = 5, dim = 256;
  SimMode update = M_JACOBI;
  char opt;

  //Default number of threads = 8
  size_t num_threads = 8;
  while ((opt = getopt(argc, argv, "hvs:r:It:g:u:")) != -1) {
    switch (opt) {
    case 'h': usage(); exit(0); break;
    case 'v': verbose = true; break;
    case 's': steps = atoi(optarg); break;
    case 'r': runs = atoi(optarg); break;
    case 'I': instrument = true; break;
    case 't': num_threads = atoi(optarg); break;
    case 'g': dim = atoi(optarg); break;
    case 'u': update = (atoi(optarg) == 0? M_JACOBI: M_REDBLACK); break;
    default:
      fprintf(stderr, "Usage: %s \n", argv[0]);
      usage();
      exit(EXIT_FAILURE);
    }
  }

  // Force the OMP to only use static scheduling
  #if defined(_OPENMP)
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(num_threads);
  #endif


  //Begin tracking the instrumentation
  track_activity(instrument); 

  //STARTUP
  start_activity(ACTIVITY_STARTUP);
  grid_t *g = new_grid(dim,dim);
  initialize_grid(g);
  finish_activity(ACTIVITY_STARTUP);

  //MAIN LOOP
  //Not parallelizable as each run is dependant on the result of the previous
  for(int i = 0; i < runs; i ++){
    if (verbose) printf("Run:\t%d\n", i);

    //Set the update 
    run_grid(g, steps, update);

    //Write_raw is used to compare the exact values of the grid
    //Whereas write_ppm only considers the concentrations of v on a scale of
    // 0-255
    //write_raw(g,i);
    write_ppm(g, i);
  }
  show_activity(instrument);
  return 0;
}

