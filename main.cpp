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

//np (number of processes is given as an arg to mpi run)
void usage(){
  printf("-h\t help\n");
  printf("-v\t verbose mode\n");
  printf("-s S\t steps per run \t\t Default: 500\n");
  printf("-I\t instrument mode \n");
  printf("-g G\t grid size \t\t Default: 256\n");
  printf("-d D\t horizantal divisions\n");
  printf("-u U\t Update function 0 for Jacobi 1 for Red Black \t Default: 0\n");
}

int main(int argc, char** argv){
  //For performance modeling
  bool instrument = false, verbose = false;

  //Default values for the grid and simulations parameters
  int steps = 500, runs = 5, gridsize = 500;
  SimMode update = M_JACOBI;
  char opt;

  //Specific MPI variables
  int process_count = 1, horizantal_divisions = 1, this_zone = 0, nzone = 0;

  //Initialize the MPI variables
  #if MPI
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);
  #endif

  bool mpi_master = this_zone == 0;

  
  while ((opt = getopt(argc, argv, "hvs:r:Ig:d:u:")) != -1) {
    switch (opt) {
    case 'h': usage(); Exit(); break;
    case 'v': verbose = true; break;
    case 's': steps = atoi(optarg);break;
    case 'r': runs = atoi(optarg);break;
    case 'I': instrument = true; break;
    case 'g': gridsize = atoi(optarg); break;
    case 'd': horizantal_divisions = atoi(optarg); break;
    default:
      fprintf(stderr, "Usage: %s \n", argv[0]);
      usage();
      Exit();
    }
  }

  //Begin tracking the instrumentation
  track_activity(instrument); 

  //STARTUP
  if (mpi_master) start_activity(ACTIVITY_STARTUP);
  state_t *s = init_zone(gridsize, gridsize, process_count, this_zone, horizantal_divisions);
  if (s == NULL) {
    fprintf(stderr, "Zone not initialized\n");
    Exit();
  }
  initialize_grid(s);
  if (mpi_master) finish_activity(ACTIVITY_STARTUP);

  //different implementation of how the grid is represented than from OMP due
  //to the lack of shared memory
  grid_t *new_g;

  //MAIN LOOP
  //Not parallelizable as each run is dependant on the result of the previous
  for(int i = 0; i < runs; i ++){
    if (verbose && mpi_master) printf("Run:\t%d\n", i);

    // the update
    new_g = run_grid(s, steps, M_REDBLACK);

    //Write_raw is used to compare the exact values of the grid
    //Whereas write_ppm only considers the concentrations of v on a scale of
    // 0-255
    //if(mpi_master) write_ppm(new_g, i);
    if(mpi_master) write_ppm(new_g, i);
  }
  if (mpi_master){
    show_activity(instrument);
  }
  Exit();
}

