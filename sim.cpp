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
      int radius = 4;
      
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
      int ind = (i + ii) * (g->ncol + 2) + j + jj;
      acc += (u? g->u : g->v)[ind] * LAPLACIAN[ii][jj];
    }
  }
  return acc;
}

void jacobi_step(grid_t *g, double *temp_u, double *temp_v){
  #pragma omp parallel
  { 
  #pragma omp for
  for(int i = 0; i < g->nrow; i += STRIDE_I){
    for(int j = 0; j < g->ncol; j += STRIDE_J){
      for(int ii = 0; ii < STRIDE_I; ii ++){
        for(int jj = 0; jj < STRIDE_J; jj ++){
          if (i + ii >=  g->nrow || j + jj >=  g->ncol ) continue;

          int ind = GINDEX(g, i + ii, j+jj);
          double Du = DU * getGrad(g, i + ii, j+jj,true);
          double Dv = DV * getGrad(g, i + ii, j+jj,false);
          double u = g->u[ind];
          double v = g->v[ind];
          u += Du - u*v*v + FEED_RATE * (1.0 - u);
          v += Dv + u*v*v - (FEED_RATE + KILL_RATE)  * v;
          u = CLAMP(u, 0.0, 1.0);
          v = CLAMP(v, 0.0, 1.0);
          
          #pragma omp atomic write
          temp_u[ind] = u;
          #pragma omp atomic write
          temp_v[ind] = v;
        }
      }
    }
  }

  }
  #pragma omp barrier

  std::swap(g->u, temp_u);
  std::swap(g->v, temp_v);
  #pragma omp barrier
}

double run_grid(grid_t *g, int steps, SimMode m) {
  double *temp_u = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  double *temp_v = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  start_activity(ACTIVITY_JSTEP);
  for(int s = 0; s < steps; s++) {
    
    jacobi_step(g, temp_u, temp_v);
  }
  finish_activity(ACTIVITY_JSTEP);
  return 1;

  free(temp_u);
  free(temp_v);
}