#include <stdio.h>
#include "math.h"
#include "stdlib.h"


void burgers1d_b(double *u, double *ub, double *u_1, double *u_1b, double D,
                 double C, int n);

void head_b (double *u, double *u_1, double C, double D, double *u_b, double *u_1_b, int n){
  burgers1d_b(u, u_b, u_1, u_1_b, D, C, n);
}

