#include "setup.h"
#include <stdlib.h>
#include <stdio.h>
//TODO make a header and also fix Makefile to include
#if MPI

state_t *divide_grid(grid_t *g, int process_count){
    int *zone_starts = (int*)malloc(sizeof(int) * process_count * 4); //Technically can be done with *2
    int index = 0;
    int row = 0;
    int col = 0;
    int proc_row = g->nrow / process_count;
    for (int p = 0; p < process_count - 1; p ++){
        zone_starts[index++] = row;
        row += proc_row;
        zone_starts[index++] = col;
        zone_starts[index++] = row;
        row += proc_row;
        zone_starts[index++] = g->ncol - 1;
    }
    zone_starts[index++] = row;
    zone_starts[index++] = row;
    zone_starts[index++] = g->nrow - 1;
    zone_starts[index++] = g->ncol - 1;
    MPI_Bcast(zone_starts, process_count * 4, MPI_INT, 0, MPI_COMM_WORLD);
    state_t *s = (state_t*)malloc(sizeof(state_t));
    s->nzones = process_count;
    s->this_zone = 0;
    s->zone_starts = zone_starts;
    s->start_row = zone_starts[0];
    s->start_col = zone_starts[1];
    s->end_row = zone_starts[2];
    s->end_col = zone_starts[3];
    return s;
}

state_t *get_divide(int process_count, int this_zone){
    state_t *s = (state_t*)malloc(sizeof(state_t));
    s->nzones = process_count;
    s->this_zone = this_zone;
    s->zone_starts = (int*)malloc(sizeof(int) * process_count * 4);
    MPI_Bcast(s->zone_starts, process_count * 4, MPI_INT, 0, MPI_COMM_WORLD);
    int index = this_zone * 4;
    s->start_row = s->zone_starts[index];
    s->start_col = s->zone_starts[index + 1];
    s->end_row = s->zone_starts[index + 2];
    s->end_col = s->zone_starts[index + 3];
    return s;
}

int get_neighbor(state_t *s, int prev_neighbor, int row, int col){
    if (prev_neighbor != -1){
        int prev_index = 4 * prev_neighbor;
        if (s->zone_starts[prev_index] < row && s->zone_starts[prev_index + 1] < col 
        &&  s->zone_starts[prev_index +2] > row && s->zone_starts[prev_index + 3] > col)
            return prev_neighbor;
    }
    for(int zid = 0; zid < s->nzones; zid ++){
        int index = zid * 4;
        if (s->zone_starts[index] < row && s->zone_starts[index + 1] < col 
        &&  s->zone_starts[index +2] > row && s->zone_starts[index + 3] > col)
            return zid;
    }
}

void init_zone(state_t *s){
    for (int zid = 0; zid < s->nzones; zid ++){
        s->transfer_count = (int *)malloc(sizeof(int) * s->nzones);
        s->import_node_list = (int **)malloc(sizeof(int *) * s->nzones);
        s->export_node_list = (int **)malloc(sizeof(int *) * s->nzones);
    }
    s->requests = (MPI_Request *)malloc(sizeof(MPI_Request) * s->nzones);
    int prev_neighbor = -1;
    for (int row = s->start_row; row < s->end_row; row ++){
        for (int col = s->start_col; col < s->end_col; col ++){
            for (int i = -1; i < 2; i ++){
                for (int j = -1; j < 2; j ++){
                    if (row + i < 0 || row + i > s->g->nrow ||
                        col + j < 0 || col + j > s->g->ncol) continue;
                    else if (row + i > s->start_row && row + i < s->end_row &&
                             col + j > s->start_col && col + j < s->end_col)
                             continue;
                    prev_neighbor = get_neighbor(s, prev_neighbor, row+i, col+j);
                    s->transfer_count[prev_neighbor]++;
                }
            }
        }
    }
    return;
    int *index = (int*)malloc(sizeof(int) * s->nzones);
    for (int zid = 0; zid < s->nzones; zid ++){
        s->import_node_list[zid] = (int *)malloc(sizeof(int) * s->transfer_count[zid] * 4);
        s->export_node_list[zid] = (int *)malloc(sizeof(int) * s->transfer_count[zid] * 4);
        index[zid]=0;
    }
     //TODO make some soft of nid so I do not keep needing to 
                 //use the start and ending row
    for (int row = s->start_row; row < s->end_row; row ++){
        for (int col = s->start_col; col < s->end_col; col ++){
            for (int i = -1; i < 2; i ++){
                for (int j = -1; j < 2; j ++){
                    if (row + i < 0 || row + i > s->g->nrow ||
                        col + j < 0 || col + j > s->g->ncol) continue;
                    else if (row + i > s->start_row && row + i < s->end_row &&
                             col + j > s->start_col && col + j < s->end_col)
                             continue;
                    prev_neighbor = get_neighbor(s, prev_neighbor, row+i, col+j);
                    s->export_node_list[prev_neighbor][index[prev_neighbor]]  = row;
                    s->export_node_list[prev_neighbor][index[prev_neighbor]+1] = col;
                    index[prev_neighbor] += 4;
                }
            }
        }
    }
}

void prepare_uv(state_t *s){
    int row, col;
    double u,v;
    for (int zid = 0; zid < s->nzones; zid ++){
        if (zid == s->this_zone) continue;
        for (int eid = 0; eid < s->transfer_count[zid]; eid ++){
            row = s->export_node_list[zid][eid*4];
            col = s->export_node_list[zid][eid*4+1];
            int index = col + row * s->g->ncol;
            u = s->g->u[index];
            v = s->g->v[index];
            s->export_node_list[zid][eid*4 + 2] = u;
            s->export_node_list[zid][eid*4 + 3] = v;
        }
    }
}

void prepare_graph(state_t *s){
    int row, col;
    double u,v;
    for (int zid = 0; zid < s->nzones; zid ++){
        if (zid == s->this_zone) continue;
        for (int eid = 0; eid < s->transfer_count[zid]; eid ++){
            row = s->import_node_list[zid][eid*4];
            col = s->import_node_list[zid][eid*4+1];
            int index = col + row * s->g->ncol;
            u = s->import_node_list[zid][eid*4 +2];
            v = s->import_node_list[zid][eid*4 +2];
            s->g->u[index] = u;
            s->g->v[index] = v;
        }
    }
}

void transfer_uv(state_t *s){
    prepare_uv(s);
    for(int zid = 0; zid < s->nzones; zid ++){
        if (zid == s->this_zone) continue;
        MPI_Isend(s->export_node_list[zid], s->transfer_count[zid],
                  MPI_INT, zid, 0, MPI_COMM_WORLD, &s->requests[zid]);
    }
    for(int zid = 0; zid < s->nzones; zid ++){
        if (zid == s->this_zone) continue;
        MPI_Recv(s->import_node_list[zid], s->transfer_count[zid],
                 MPI_INT, zid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (int zid = 0; zid < s->nzones; zid ++){
        if (zid == s->this_zone) continue;
        MPI_Wait(&s->requests[zid], MPI_STATUS_IGNORE);
    }
    prepare_graph(s);
}
#endif