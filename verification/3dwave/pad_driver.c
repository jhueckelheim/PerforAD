#include <stdio.h>
#include "math.h"
#include "stdlib.h"


void wave3d_perf_b(double*** u, double*** u_b, double*** c, double*** u_1, double*** u_1_b, double*** u_2, double*** u_2_b, double D, int n) ;

void head_b(double*** u, double*** u_b, double*** c, double*** u_1, double*** u_1_b, double*** u_2, double*** u_2_b, double D, int n) {
  wave3d_perf_b(u, u_b, c, u_1, u_1_b, u_2, u_2_b, D, n);

}
