#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <utility>
#include <algorithm>
#include "sim.h"
#include "setup.h"

grid_t *new_grid(int nrow, int ncol) {
  grid_t *g = (grid_t*) malloc(sizeof(grid_t));
  g->nrow = nrow;
  g->ncol = ncol;
  g->u = (double*) calloc((nrow+2) * (ncol+2), sizeof(double));
  g->v = (double*) calloc((nrow+2) * (ncol+2), sizeof(double));

  return g;
}

void free_grid(grid_t *g) {
  free(g->u);
  free(g->v);
  free(g);
}

void initialize_grid(state_t *s, InitMode m) {
  grid_t *g = s->g;
  for(int i = -1; i <= g->nrow; i++) {
    for(int j = -1; j <= g->ncol; j++) {
      int ind = GINDEX(g,i,j);

      int real_i = i + s->start_row;
      int real_j = j + s->start_col;
      //border
      int centery = s->nrow/2;
      int centerx = s->ncol/2;
      int radius = centerx / 4;
      
      if ((real_i-centery) * (real_i-centery) + (real_j-centerx) * (real_j-centerx) <= radius * radius){
        g->u[ind] = 1.0/4;
        g->v[ind] = 1.0/2;
      } else {
        g->u[ind] = 1.0;
        g->v[ind] = 0;
      }
    }
  }
}

double getGrad(grid_t *g, int i, int j, bool u) {
  double acc = 0;
  for(int ii = 0; ii < K_SIZE; ii++) {
    for(int jj = 0; jj < K_SIZE; jj++) {
      if(LAPLACIAN[ii][jj] == 0) continue;
      int ind = (i + ii) * (g->ncol + 2) + j + jj;
      acc += (u? g->u : g->v)[ind] * LAPLACIAN[ii][jj];
    }
  }
  return acc;
}

void jacobi_step(state_t *s){
  if (s->this_zone == 0) start_activity(ACTIVITY_SSTEP);

  grid_t *g = s->g;
  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j ++){
      int ind = GINDEX(g, i, j);
      double Du = DU * getGrad(g, i, j,true);
      double Dv = DV * getGrad(g, i, j,false);
      double u = g->u[ind];
      double v = g->v[ind];
      u += Du - u*v*v + FEED_RATE * (1.0 - u);
      v += Dv + u*v*v - (FEED_RATE + KILL_RATE)  * v;
      u = CLAMP(u, 0.0, 1.0);
      v = CLAMP(v, 0.0, 1.0);
      
      (g->temp_u)[ind] = u;
      (g->temp_v)[ind] = v;
    }
  }
  std::swap(g->u, g->temp_u);
  std::swap(g->v, g->temp_v);

  if (s->this_zone == 0) finish_activity(ACTIVITY_SSTEP);
}

void red_black_step(state_t *s, int parity){
  if (s->this_zone == 0) start_activity(ACTIVITY_SSTEP);
 
  grid_t *g = s->g;
  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j ++){
      if((s->start_row + i + s->start_col + j) % 2 == parity) continue;
      int ind = GINDEX(g, i, j);
      double Du = DU * getGrad(g, i, j,true);
      double Dv = DV * getGrad(g, i, j,false);
      double u = g->u[ind];
      double v = g->v[ind];
      u += Du - u*v*v + FEED_RATE * (1.0 - u);
      v += Dv + u*v*v - (FEED_RATE + KILL_RATE)  * v;
      u = CLAMP(u, 0.0, 1.0);
      v = CLAMP(v, 0.0, 1.0);
      
      g->u[ind] = u;
      g->v[ind] = v;
    }
  } 
  
  if (s->this_zone == 0) finish_activity(ACTIVITY_SSTEP);
}

grid_t *run_grid(state_t *state, int steps, SimMode m) {
  grid_t *g = state->g;

  if(m == M_JACOBI) {
    g->temp_u = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
    g->temp_v = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  }

  for(int s = 0; s < steps; s++) {    
    if(m == M_JACOBI) {
      jacobi_step(state);
      #if MPI
        begin_exchange_uv(state);
        finish_exchange_uv(state);
      #endif
    } else if (m == M_REDBLACK) {
      red_black_step(state, 0);
      #if MPI
        begin_exchange_uv(state);
        finish_exchange_uv(state);
      #endif
      red_black_step(state, 1);
      #if MPI
        begin_exchange_uv(state);
        finish_exchange_uv(state);
      #endif
    }
  }
  //send all to master
  grid_t *new_g;
  #if MPI
  if (state->this_zone == 0) start_activity(ACTIVITY_GLOBAL_COMM);
  if(state->this_zone == 0){
    new_g = gather_uv(state);
  } else {
    send_uv(state);
  }
  if (state->this_zone == 0) finish_activity(ACTIVITY_GLOBAL_COMM);
  #else
  new_g = state->g;
  #endif
  return new_g;
}