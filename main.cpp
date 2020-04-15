#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string>
#include "sim.h"

void print_u(grid_t *g, int iter){
  char buf[50];
  sprintf(buf, "out%d.ppm", iter);
  FILE *fp = fopen(buf, "wb");

  if (!fp) {
      fprintf(stderr, "Error: could not open file for write\n");
      exit(1);
  }
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", g->ncol, g->nrow);
  fprintf(fp, "255\n");
  for (int j=g->nrow-1; j>=0; j--) {
    for (int i=0; i< g->ncol; i++) {

      double dub = g->v[GINDEX(g,j,i)] * 3;

      char val[3];
      val[0] = static_cast<char>(255. * CLAMP(dub, 0., 1.));
      val[1] = static_cast<char>(255. * CLAMP(dub, 0., 1.));
      val[2] = static_cast<char>(255. * CLAMP(dub, 0., 1.));

      fputc(val[0], fp);
      fputc(val[1], fp);
      fputc(val[2], fp);
    }
  }

    fclose(fp);

}

int main(int argc, char** argv)
{
  //args: Time steps

  grid_t *g = new_grid(500,500);
  initialize_grid(g);
  for(int i = 0; i < 20; i ++){
    run_grid(g, 500);
    print_u(g, i);
  }
  return 0;
}

