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

void initialize_grid(grid_t *g, InitMode m) {
  for(int i = -1; i <= g->nrow; i++) {
    for(int j = -1; j <= g->ncol; j++) {
      int ind = GINDEX(g,i,j);
      //border
      int centery = g->nrow/2;
      int centerx = g->ncol/2;
      int radius = centerx / 4;
      
      if ((i-centery) * (i-centery) + (j-centerx) * (j-centerx) <= radius * radius){
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
  grid_t *g = s->g;
  for(int i = s->start_row; i < s->end_row; i ++){
    for(int j = s->start_col; j < s->end_col; j ++){
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
}

void red_black_step(state_t *s){
  grid_t *g = s->g;
  for(int parity = 0; parity < 2; parity++){
    for(int i = s->start_row; i < s->end_row; i ++){
      for(int j = s->start_col; j < s->end_col; j ++){
        if((i+j) % 2 == parity) continue;
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
  }
}

double run_grid(state_t *state, int steps, SimMode m) {
  grid_t *g = state->g;

  if(m == M_JACOBI) {
    g->temp_u = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
    g->temp_v = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  }

  for(int s = 0; s < steps; s++) {
    if (state->this_zone == 0) start_activity(ACTIVITY_JSTEP);
    
    if(m == M_JACOBI) {
      jacobi_step(state);
    } else if (m == M_REDBLACK) {
      red_black_step(state);
    }
    
    
    if (state->this_zone == 0) finish_activity(ACTIVITY_JSTEP);
    //exchange here
    #if MPI
      if (state->this_zone == 0) start_activity(ACTIVITY_LOCAL_COMM);
      exchange_uv(state);
      if (state->this_zone == 0) finish_activity(ACTIVITY_LOCAL_COMM);
    #endif 
  }
  //send all to master
  #if MPI
  if (state->this_zone == 0) start_activity(ACTIVITY_GLOBAL_COMM);
  if(state->this_zone == 0){
    gather_uv(state);
  } else {
    send_uv(state);
  }
  if (state->this_zone == 0) finish_activity(ACTIVITY_GLOBAL_COMM);
  #endif
  return 0;
}