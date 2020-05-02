#ifndef SIM_H
#define SIM_H

#include "instrument.h"
#ifndef MPI
#define MPI 0
#endif

#if MPI
  #include <mpi.h>
#endif

#define STRIDE_I 1
#define STRIDE_J 1

#define GINDEX(g,r,c) (((r)+1)*((g)->ncol+2)+((c)+1))

#define CLAMP(x, minimum, maximum) std::max(minimum, std::min(x, maximum))

enum InitMode {UNIFORM};
enum SimMode {M_JACOBI, M_REDBLACK};
enum Direction {NORTH = 0, SOUTH, EAST, WEST};

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

  double *temp_u; //padded +2
  double *temp_v; //padded +2
} grid_t;

typedef struct{
  grid_t *g;
  int nzones;
  int h_divs;
  /*** Info on this specific ZONE ***/
  int this_zone;

  //Information pertaining to the partitioning of the zones.
  int start_row;
  int start_col;
  int end_row;
  int end_col;

  //The number of neighbors and which nodes we need to send and recieve from them
  int *neighbors;
  //Should only count for the ones this zone actually has neighbors with.
  //there will be a difference between the neghibors zid and their location
  //in the import and export lists 
  //TODO make this less confusing
  double **export_node_list;
  double **import_node_list;

  //Keep trake of the MPI requests after the sending
  #if MPI
    MPI_Request *send_requests;
    MPI_Request *recv_requests;
  #endif
} state_t;

grid_t *new_grid(int nrow, int ncol);
void free_grid(grid_t *g);
void initialize_grid(grid_t *g, InitMode m = UNIFORM);

double run_grid(state_t *s, int steps, SimMode m = M_REDBLACK);

#endif