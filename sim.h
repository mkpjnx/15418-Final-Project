#include "instrument.h"

#define STRIDE_I 3
#define STRIDE_J 32

enum InitMode {UNIFORM};
enum SimMode {JACOBI};

static const int K_SIZE = 3;

static const double FEED_RATE = 0.078;
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

#define GINDEX(g,r,c) (((r)+1)*((g)->ncol+2)+((c)+1))

#define CLAMP(x, minimum, maximum) std::max(minimum, std::min(x, maximum))

grid_t *new_grid(int nrow, int ncol);
void free_grid(grid_t *g);
void initialize_grid(grid_t *g, InitMode m = UNIFORM);

double run_grid(grid_t *g, int steps, SimMode m = JACOBI);

