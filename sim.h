#ifndef SIM_H
#define SIM_H

#include "instrument.h"
#ifndef MPI
#define MPI 0
#endif

#if MPI
  #include <mpi.h>
#endif

#define STRIDE_I 3
#define STRIDE_J 32

enum InitMode {UNIFORM};
enum SimMode {JACOBI};

static const int K_SIZE = 3;

static const double FEED_RATE = 0.037;
static const double KILL_RATE = 0.0612;
static const double DU = 0.209;
static const double DV = 0.105;


static const double LAPLACIAN[K_SIZE][K_SIZE] =
  {{0,1,0},{1,-4,1},{0,1,0}};

typedef struct {
  int nrow;
  int ncol;
  double *u; //padded +2
  double *v; //padded +2
} grid_t;

typedef struct{
  grid_t *g;
  int nzones;
  /*** Info on this specific ZONE ***/
  int this_zone;

  //Information pertaining to the partitioning of the zones.
  int* zone_starts;
  int start_row;
  int start_col;
  int end_row;
  int end_col;

  //The number of neighbors and which nodes we need to send and recieve from them
  int *transfer_count;
  //Should only count for the ones this zone actually has neighbors with.
  //there will be a difference between the neghibors zid and their location
  //in the import and export lists 
  //TODO make this less confusing
  int **export_node_list;
  int **import_node_list;

  //Keep trake of the MPI requests after the sending
  #if MPI
    MPI_Request *requests;
  #endif;
} state_t;
#define GINDEX(g,r,c) (((r)+1)*((g)->ncol+2)+((c)+1))

#define CLAMP(x, minimum, maximum) std::max(minimum, std::min(x, maximum))

grid_t *new_grid(int nrow, int ncol);
void free_grid(grid_t *g);
void initialize_grid(grid_t *g, InitMode m = UNIFORM);

double run_grid(state_t *s, int steps, SimMode m = JACOBI);

#endif