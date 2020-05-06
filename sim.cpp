#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <utility>
#include <algorithm>
#include "sim.h"

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

double inline getGrad(grid_t *g, int i, int j, bool u) {
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

void inline calc_single(grid_t *g, int i, int j, bool in_place = true) {
    int ind = GINDEX(g, i, j);
    double Du = DU * getGrad(g, i, j,true);
    double Dv = DV * getGrad(g, i, j,false);
    double u = g->u[ind];
    double v = g->v[ind];
    u += Du - u*v*v + FEED_RATE * (1.0 - u);
    v += Dv + u*v*v - (FEED_RATE + KILL_RATE)  * v;
    u = CLAMP(u, 0.0, 1.0);
    v = CLAMP(v, 0.0, 1.0);
    
    (in_place ? g->u : g->temp_u)[ind] = u;
    (in_place ? g->v : g->temp_v)[ind] = v;
}


void jacobi_step(grid_t *g){
  #pragma omp parallel
  { 
  #pragma omp for
  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j ++){
      calc_single(g, i ,j, false);
    }
  }
  }
  std::swap(g->u, g->temp_u);
  std::swap(g->v, g->temp_v);
}

void red_black_step(grid_t *g){
  #pragma omp parallel
  { 
  #pragma omp for
  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j ++){
      if((i+j)%2 == 0) continue;
      calc_single(g,i,j);
      }
    }

  #pragma omp for
  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j ++){
      if((i+j)%2 == 1) continue;
      calc_single(g,i,j);
    }
  }

  }
}

double run_grid(grid_t *g, int steps, SimMode m) {
  g->temp_u = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  g->temp_v = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  start_activity(ACTIVITY_JSTEP);
  for(int s = 0; s < steps; s++) {
    if (m == M_JACOBI) {
      jacobi_step(g);
    } else if (m == M_REDBLACK) {
      red_black_step(g);
    }
  }
  finish_activity(ACTIVITY_JSTEP);
  return 1;
}