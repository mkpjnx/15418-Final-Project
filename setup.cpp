#include "setup.h"
#include "sim.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <utility>
#include <algorithm>

//helper routine for calculating length of 
static inline int get_length(state_t *s, int dir) {
  return (dir == NORTH || dir == SOUTH) ? s->end_col - s->start_col : s->end_row - s->start_row;
}

/*

*|   N O R T H    |*
--------------------
W| this zone ...  |E
E| ....           |A
T| ....           |S
S| ....           |T
--------------------
*|   S O U T H    |*

*/

state_t *init_zone(grid_t *g, int process_count, int this_zone, int h_divs) {
  if (process_count % h_divs) return NULL; //unable to divide into grids

  //round up padding
  int z_width = (g->ncol + h_divs - 1) / h_divs;
  int v_divs = process_count / h_divs;
  int z_height = (g->nrow + v_divs - 1) / v_divs;

  state_t *s = (state_t*)malloc(sizeof(state_t));
  s->g = g;
  s->h_divs = h_divs;
  s->this_zone = this_zone;
  s->nzones = process_count;

  //calculate rectangle
  s->start_row = z_height * (this_zone / h_divs);
  s->start_col = z_width * (this_zone % h_divs);
  s->end_col = std::min(s->start_col + z_width, g->ncol);
  s->end_row = std::min(s->start_row + z_height, g->nrow);

  //Check if there are neighbor zones
  s->neighbors = (int*)malloc(sizeof(int) * 4);
  s->neighbors[NORTH] = this_zone / h_divs ? this_zone - h_divs : -1;
  s->neighbors[SOUTH] = this_zone / h_divs < v_divs - 1 ? this_zone + h_divs : -1;
  s->neighbors[EAST] = this_zone % h_divs < h_divs - 1 ? this_zone + 1 : -1;
  s->neighbors[WEST] = this_zone % h_divs ? this_zone - 1 : -1;
  
  s->export_node_list = (double**)calloc(sizeof(double*), 4);
  s->import_node_list = (double**)calloc(sizeof(double*), 4);

  //allocate buffers
  for(int dir = 0; dir < 4; dir ++) {
    if(s->neighbors[dir] == -1) continue;
    int length = get_length(s, dir);
    s->export_node_list[dir] = (double*)malloc(sizeof(double) * length * 4);
    s->import_node_list[dir] = s->export_node_list[dir] + 2 * length;
  }
  
  #if MPI
  s->requests = (MPI_Request *)malloc(sizeof(MPI_Request) * 4);
  #endif
  return s;
}

#if MPI

void exchange_uv(state_t *s){
  grid_t *g = s->g;

  //export step
  for(int dir = 0; dir < 4; dir ++) {
    int zone = s->neighbors[dir];
    if(zone == -1) continue;

    int length = get_length(s, dir);
    //load into export buffer
    for(int k = 0; k < length; k++) {
      int grid_ind;
      switch(dir){
        case NORTH: grid_ind = GINDEX(g, s->start_row, s->start_col + k); break;
        case SOUTH: grid_ind = GINDEX(g, s->end_row - 1, s->start_col + k); break;
        case EAST:  grid_ind = GINDEX(g, s->start_row + k, s->end_col - 1); break;
        case WEST:  grid_ind = GINDEX(g, s->start_row + k, s->start_col); break;
      }

      s->export_node_list[dir][k] = g->u[grid_ind];
      s->export_node_list[dir][k + length] = g->v[grid_ind];
      //if(s->this_zone == 0) printf("dir %d\t grid_ind %d\n", dir, grid_ind);
    }
    //async send
    MPI_Isend(s->export_node_list[dir], 2*length, MPI_DOUBLE, zone, 0,
            MPI_COMM_WORLD, &s->requests[dir]);
  }

  //recv
  for(int dir = 0; dir < 4; dir ++) {
    int zone = s->neighbors[dir];
    if(zone == -1) continue;
    int length = get_length(s, dir);
    MPI_Recv(s->import_node_list[dir], 2*length, MPI_DOUBLE, zone, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //wait
  for(int dir = 0; dir < 4; dir ++) {
    int zone = s->neighbors[dir];
    if(zone == -1) continue;
    MPI_Wait(&s->requests[dir], MPI_STATUS_IGNORE);
  }
  
  //load into g
  for(int dir = 0; dir < 4; dir ++) {
    int zone = s->neighbors[dir];
    if(zone == -1) continue;
    int length = get_length(s, dir);
    //load import buffer
    for(int k = 0; k < length; k++) {
      int grid_ind;
      switch(dir){
        case NORTH: grid_ind = GINDEX(g, s->start_row - 1, s->start_col + k); break;
        case SOUTH: grid_ind = GINDEX(g, s->end_row, s->start_col + k); break;
        case EAST:  grid_ind = GINDEX(g, s->start_row + k, s->end_col); break;
        case WEST:  grid_ind = GINDEX(g, s->start_row + k, s->start_col - 1); break;
      }
      g->u[grid_ind] = s->import_node_list[dir][k];
      g->v[grid_ind] = s->import_node_list[dir][k + length];
    }
  }
}

//gather after all runs are completed
void gather_uv(state_t *s){
  int h_divs = s->h_divs;
  grid_t *g = s->g;
  int z_width = (g->ncol + h_divs - 1) / h_divs;
  int v_divs = s->nzones / h_divs;
  int z_height = (g->nrow + v_divs - 1) / v_divs;
  double *buf = (double*)malloc(sizeof(double) * z_width * z_height * 2);

  for(int z = 1; z < s->nzones; z++) {
    int start_row = z_height * (z / h_divs);
    int start_col = z_width * (z % h_divs);
    int end_col = std::min(start_col + z_width, g->ncol);
    int end_row = std::min(start_row + z_height, g->nrow);

    int length = (end_col - start_col) * (end_row - start_row);
    MPI_Recv(buf, 2*length, MPI_DOUBLE, z, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int buf_ind = 0;
    for(int i = start_row; i < end_row; i++) {
      for(int j = start_col; j < end_col; j++){
        g->u[GINDEX(g,i,j)] = buf[buf_ind];
        g->v[GINDEX(g,i,j)] = buf[buf_ind + length];
        buf_ind ++;
      }
    }
  }
}

//send to master
void send_uv(state_t *s){
  int length = (s->end_col - s->start_col) * (s->end_row - s->start_row);
  double *buf = (double*)malloc(sizeof(double) * length * 2);
  int buf_ind = 0;
  MPI_Request request;
  for(int i = s->start_row; i < s->end_row; i++) {
    for(int j = s->start_col; j < s->end_col; j++){
      buf[buf_ind] = s->g->u[GINDEX(s->g,i,j)];
      buf[buf_ind + length] = s->g->v[GINDEX(s->g,i,j)];
      buf_ind ++;
    }
  }
  MPI_Isend(buf, 2*length, MPI_DOUBLE, 0, 0,
            MPI_COMM_WORLD, &request);
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  free(buf);
}
#endif