#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <utility>
#include <algorithm>
#include "sim.h"
#include "setup.h"

//Allocates enough space for local grid
grid_t *new_grid(int nrow, int ncol) {
  grid_t *g = (grid_t*) malloc(sizeof(grid_t));
  g->nrow = nrow;
  g->ncol = ncol;
  g->u = (double*) calloc((nrow+2) * (ncol+2), sizeof(double));
  g->v = (double*) calloc((nrow+2) * (ncol+2), sizeof(double));

  return g;
}

//Frees grid once it is no longer in use
void free_grid(grid_t *g) {
  free(g->u);
  free(g->v);
  free(g);
}

//Initializes the grid to have a circle of u = .25 and v = .5 centered around the 
//middle
//Sets the rest of the grid to have u = 1 and v = 0
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

//Writes the values of u and v to an output file
//Only used for checking correctess
void write_raw(grid_t *g, int iter){
  char buf[50];
  sprintf(buf, "out/out%d.txt", iter);
  FILE *fp = fopen(buf, "wb");

  if (!fp) {
      fprintf(stderr, "Error: could not open file for write\n");
      exit(1);
  }
  for (int j=g->nrow-1; j>=0; j--) {
    for (int i=0; i< g->ncol; i++) {

      double dub = g->v[GINDEX(g,j,i)];
      fprintf(fp, "%d,%d: %a\n", j, i, dub);
    }
  }
    fclose(fp);
}

//Writes the value of V clamped to 0-255 into a ppm file.
void write_ppm(grid_t *g, int iter){
  char buf[50];
  #if MPI
    sprintf(buf, "out/MPI%d.ppm", iter);
  #else
    sprintf(buf, "out/out%d.ppm", iter);
  #endif
  FILE *fp = fopen(buf, "wb");

  if (!fp) {
      fprintf(stderr, "Error: could not open file for write\n");
      exit(1);
  }

  //print header for the PPM file. has ncols and nrows
  fprintf(fp, "P6\n%d %d\n255\n", g->ncol, g->nrow);
  
  //write out the value of the grid
  for (int j=g->nrow-1; j>=0; j--) {
    for (int i=0; i< g->ncol; i++) {
      double dub = g->v[GINDEX(g,j,i)] * 3;
      char value = static_cast<char>(255. * CLAMP(dub, 0., 1.));
      char val[3] = {value, value, value};
      for (int rgb = 0; rgb < 3; rgb ++){
        fputc(val[rgb], fp);
      }
    }
  }
    fclose(fp);
}

//Loops over the Laplacian defined in sim.h to get the surrounding 
//concentrations of u and v
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

//Used to calculate the change in each grid point over a single time step.
//GINDEX is used to make sure that edge cases are avoided
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

//Jacobi update - create a dummy grid so store updates and then load that an the
//new grid
void jacobi_step(state_t *s){

  #if MPI
    begin_exchange_uv(s);
  #endif

  if (s->this_zone == 0) start_activity(ACTIVITY_SSTEP);

  grid_t *g = s->g;
  for(int i = 1; i < g->nrow - 1; i ++){
    for(int j = 1; j < g->ncol - 1; j ++){
      calc_single(g, i ,j, false);
    }
  }
  if (s->this_zone == 0) finish_activity(ACTIVITY_SSTEP);


  #if MPI
    finish_exchange_uv(s);
  #endif

  //borders
  if (s->this_zone == 0) start_activity(ACTIVITY_SSTEP);
  for (int i = 0; i < g->nrow; i += g->nrow-1) {
    for(int j = 1; j < g->ncol - 1; j ++){
      calc_single(g, i ,j, false);
    }
  }

  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j += g->ncol - 1){
      calc_single(g, i ,j, false);
    }
  }
  if (s->this_zone == 0) finish_activity(ACTIVITY_SSTEP);
  std::swap(g->u, g->temp_u);
  std::swap(g->v, g->temp_v);

}

//Red Black update - first update the red nodes (even graph index) then 
//update the black nodes (odd graph index)
void red_black_step(state_t *s, int parity){
  #if MPI
    begin_exchange_uv(s);
  #endif

  if (s->this_zone == 0) start_activity(ACTIVITY_SSTEP);
 
  grid_t *g = s->g;
  for(int i = 1; i < g->nrow - 1; i ++){
    for(int j = 1; j < g->ncol - 1; j ++){
      if((s->start_row + i + s->start_col + j) % 2 == parity) continue;
      calc_single(g,i,j);
    }
  } 
  
  if (s->this_zone == 0) finish_activity(ACTIVITY_SSTEP);

  #if MPI
    finish_exchange_uv(s);
  #endif

  if (s->this_zone == 0) start_activity(ACTIVITY_SSTEP);
  for (int i = 0; i < g->nrow; i += g->nrow-1) {
    for(int j = 1; j < g->ncol - 1; j ++){
      if((s->start_row + i + s->start_col + j) % 2 == parity) continue;
      calc_single(g, i ,j);
    }
  }

  for(int i = 0; i < g->nrow; i ++){
    for(int j = 0; j < g->ncol; j += g->ncol - 1){
      if((s->start_row + i + s->start_col + j) % 2 == parity) continue;
      calc_single(g, i ,j);
    }
  }
  
  if (s->this_zone == 0) finish_activity(ACTIVITY_SSTEP);
}

//Function to run s update steps

grid_t *run_grid(state_t *state, int steps, SimMode m) {
  grid_t *g = state->g;

  if(m == M_JACOBI) {
    g->temp_u = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
    g->temp_v = (double*) calloc((g->nrow+2) * (g->ncol+2), sizeof(double));
  }

  for(int s = 0; s < steps; s++) {    
    if(m == M_JACOBI) {
      jacobi_step(state);
    } else if (m == M_REDBLACK) {
      //red
      red_black_step(state, 0);
      //Comunicate all needed rows/cols
      #if MPI
        begin_exchange_uv(state);
        finish_exchange_uv(state);
      #endif
      //black
      red_black_step(state, 1);
      #if MPI
      //comunicate all needed rows/cols
        begin_exchange_uv(state);
        finish_exchange_uv(state);
      #endif
    }
  }
  //initialize return value
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