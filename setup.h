#include "sim.h"

#if MPI
state_t *divide_grid(grid_t *g, int process_count);
state_t *get_divide(int process_count, int this_zone);
int get_neighbor(state_t *s, int prev_neighbor, int row, int col);
void init_zone(state_t *s);
void send_uv(state_t *s);
void recv_uv(state_t *s);
#endif
