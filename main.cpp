#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string>
#include "sim.h"
#define CLAMP(x, minimum, maximum) std::max(minimum, std::min(x, maximum))

void print_u(grid_t *g, double max_u){
  FILE *fp = fopen("out.ppm", "wb");

  if (!fp) {
      fprintf(stderr, "Error: could not open file for write\n");
      exit(1);
  }
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", g->ncol, g->nrow);
  fprintf(fp, "255\n");
  for (int j=g->nrow-1; j>=0; j--) {
    for (int i=0; i< g->ncol; i++) {

      double dub = g->u[GINDEX(g,j,i)] / max_u;

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
  grid_t *g = new_grid(1000,1000);
  initialize_grid(g);
  double max_u = run_grid(g, 1000);

  print_u(g, max_u);
  return 0;
}

