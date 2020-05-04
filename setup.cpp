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

static inline void get_zone_start_end(int total, int divs, int ind, int *start, int *end){
  int base = total/divs;
  int big = total % divs;
  *start = (ind > big) ? base * ind + big : (base+1) * ind;
  *end = *start + base + ((ind < big) ? 1 : 0);
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

state_t *init_zone(int nrow, int ncol, int process_count, int this_zone, int h_divs) {
  if (process_count % h_divs) return NULL; //unable to divide into grids

  //round up padding
  int z_width = ncol / h_divs;
  int v_divs = process_count / h_divs;
  int z_height = nrow / v_divs;

  state_t *s = (state_t*)malloc(sizeof(state_t));
  s->nrow = nrow;
  s->ncol = ncol;

  s->h_divs = h_divs;
  s->this_zone = this_zone;
  s->nzones = process_count;

  //calculate rectangle
  get_zone_start_end(nrow, v_divs, this_zone / h_divs,
    &(s->start_row), &(s->end_row));
  
  get_zone_start_end(ncol, h_divs, this_zone % h_divs,
    &(s->start_col), &(s->end_col));

  //local grid
  s->g = new_grid(s->end_row - s->start_row, s->end_col - s->start_col);

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
  s->send_requests = (MPI_Request *)malloc(sizeof(MPI_Request) * 4);
  s->recv_requests = (MPI_Request *)malloc(sizeof(MPI_Request) * 4);
  #endif
  return s;
}

#if MPI


void begin_exchange_uv(state_t *s){
  if (s->this_zone == 0) start_activity(ACTIVITY_LOCAL_COMM);

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
        case NORTH: grid_ind = GINDEX(g, 0, k); break;
        case SOUTH: grid_ind = GINDEX(g, g->nrow - 1, k); break;
        case EAST:  grid_ind = GINDEX(g, k, g->ncol - 1); break;
        case WEST:  grid_ind = GINDEX(g, k, 0); break;
      }

      s->export_node_list[dir][k] = g->u[grid_ind];
      s->export_node_list[dir][k + length] = g->v[grid_ind];
    }
    //async send
    MPI_Isend(s->export_node_list[dir], 2*length, MPI_DOUBLE, zone, 0,
            MPI_COMM_WORLD, &s->send_requests[dir]);
  }

  //async recv
  for(int dir = 0; dir < 4; dir ++) {
    int zone = s->neighbors[dir];
    if(zone == -1) continue;
    int length = get_length(s, dir);
    MPI_Irecv(s->import_node_list[dir], 2*length, MPI_DOUBLE, zone, 0,
            MPI_COMM_WORLD, &s->recv_requests[dir]);
  }

  if (s->this_zone == 0) finish_activity(ACTIVITY_LOCAL_COMM);
}


void finish_exchange_uv(state_t *s){
  if (s->this_zone == 0) start_activity(ACTIVITY_LOCAL_COMM);
  
  grid_t *g = s->g;

  //wait
  for(int dir = 0; dir < 4; dir ++) {
    int zone = s->neighbors[dir];
    if(zone == -1) continue;
    MPI_Wait(&s->send_requests[dir], MPI_STATUS_IGNORE);
    MPI_Wait(&s->recv_requests[dir], MPI_STATUS_IGNORE);
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
        case NORTH: grid_ind = GINDEX(g, -1, k); break;
        case SOUTH: grid_ind = GINDEX(g, g->nrow, k); break;
        case EAST:  grid_ind = GINDEX(g, k, g->ncol); break;
        case WEST:  grid_ind = GINDEX(g, k, -1); break;
      }
      g->u[grid_ind] = s->import_node_list[dir][k];
      g->v[grid_ind] = s->import_node_list[dir][k + length];
    }
  }

  if (s->this_zone == 0) finish_activity(ACTIVITY_LOCAL_COMM);
}

//gather after all runs are completed
void gather_uv(state_t *s){
  int h_divs = s->h_divs;
  grid_t *g = new_grid(s->nrow, s->ncol);
  int z_width = s->ncol / h_divs;
  int v_divs = s->nzones / h_divs;
  int z_height = s->nrow / v_divs;
  double *buf = (double*)malloc(sizeof(double) * (z_width+1)*(z_height+1)*2);

  //from self
  for(int i = 0; i < s->g->nrow; i++) {
    for(int j = 0; j < s->g->ncol; j++){
      g->u[GINDEX(g,i + s->start_row,j + s->start_col)] = s->g->u[GINDEX(s->g,i,j)];
      g->v[GINDEX(g,i + s->start_row,j + s->start_col)] = s->g->v[GINDEX(s->g,i,j)];
    }
  }

  //from others
  for(int z = 1; z < s->nzones; z++) {
    int start_row, start_col, end_col, end_row;
    //calculate rectangle
    get_zone_start_end(s->nrow, v_divs, z / h_divs,
      &(start_row), &(end_row));
    
    get_zone_start_end(s->ncol, h_divs, z % h_divs,
      &(start_col), &(end_col));


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
  s->g = g;
}

//send to master
void send_uv(state_t *s){
  int length = (s->end_col - s->start_col) * (s->end_row - s->start_row);
  double *buf = (double*)malloc(sizeof(double) * length * 2);
  int buf_ind = 0;
  MPI_Request request;
  for(int i = 0; i < s->g->nrow; i++) {
    for(int j = 0; j < s->g->ncol; j++){
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