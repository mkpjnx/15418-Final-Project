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
  for(int i = 0; i < g->nrow + 2; i++) {
    for(int j = 0; j < g->ncol + 2; j++) {
      int ind = i * (g->ncol + 2) + j;
      //border
      if(!i || i == g->nrow + 1 || !j || j == g->ncol + 1) {
        g->u[ind] = 1.;
      } else {
         g->u[ind] = .5;
         g->v[ind] = .25;
      }
    }
  }
}

double getGrad(grid_t *g, int i, int j, bool u) {
  double acc = 0;
  for(int ii = 0; ii < K_SIZE; ii++) {
    for(int jj = 0; jj < K_SIZE; jj++) {
      acc += (u? g->u : g->v)[GINDEX(g, i+ii-1, j+jj-1)] * LAPLACIAN[ii][jj];
    }
  }

  return acc;
}

static inline double clamp(double n) {
  return std::max(0.0, std::min(n, 1.0));
}

void jacobi_step(grid_t *g, double *temp_u, double *temp_v){
  for(int i = 0; i < g->nrow; i++){
      for(int j = 0; j < g->ncol; j++){
        int ind = GINDEX(g, i, j);
        double Du = getGrad(g,i,j,true);
        double Dv = getGrad(g,i,j,false);
        double u = g->u[ind];
        double v = g->v[ind];

        temp_u[ind] = clamp(Du - u*v*v + FEED_RATE * (1 - u));
        temp_v[ind] = clamp(Dv + u*v*v - (FEED_RATE + KILL_RATE)  * v);
        
    }
  }
  std::swap(g->u, temp_u);
  std::swap(g->v, temp_v);
}

double run_grid(grid_t *g, int steps, SimMode m) {
  double *temp_u = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  double *temp_v = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  
  for(int s = 0; s < steps; s++) {
    printf("step %d: %f\n", s, g->u[GINDEX(g,0,0)]);
    jacobi_step(g, temp_u, temp_v);
  }
  double max_u = 0;
  for(int i = 0; i < g->nrow; i++){
    for(int j = 0; j < g->ncol; j++){
      int ind = GINDEX(g, i, j);
      max_u = std::max(g->u[ind], max_u);
    }
  }
  return max_u;
}