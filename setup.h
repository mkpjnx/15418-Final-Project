#include "sim.h"


state_t *init_zone(grid_t *g, int process_count, int this_zone, int h_divs);
#if MPI
void begin_exchange_uv(state_t *s);

void finish_exchange_uv(state_t *s);

void gather_uv(state_t *s);

void send_uv(state_t *s);
#endif
