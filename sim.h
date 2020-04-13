
enum InitMode {UNIFORM};
enum SimMode {JACOBI};

static const int K_SIZE = 3;

static const double FEED_RATE = 0.01;
static const double KILL_RATE = 0.047;

static const double LAPLACIAN[K_SIZE][K_SIZE] =
  {{0,.25,0},{0.25,-1,0.25},{0,.25,0}};

typedef struct {
  int nrow;
  int ncol;
  double *u; //padded +2
  double *v; //padded +2
} grid_t;

#define GINDEX(g,r,c) (((r)+1)*((g)->ncol+2)+((c)+1))

grid_t *new_grid(int nrow, int ncol);
void free_grid(grid_t *g);
void initialize_grid(grid_t *g, InitMode m = UNIFORM);

double run_grid(grid_t *g, int steps, SimMode m = JACOBI);

