#include <stdio.h>
#include "math.h"
#include "stdlib.h"


void burgers1d_perf_b(double* u, double* u_b, double* u_1, double* u_1_b, double D, double C, int n);
/*
#ifndef _CIVL
void head_b (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b, int n);
int main()
{
  double outv[N], inv[N], vel[N];
  double outv_b[N], inv_b[N], vel_b[N];
  int i,j;
  for (i = 0; i < N; i++){
    outv[i] = 0.0;
    inv[i] = (i+1)*(i+1) * 0.5;
    vel[i] = (i+1)*(i+1) * 0.5;
    outv_b[i] = 1.0;
    inv_b[i] = 0.0;
    vel_b[i] = 0.0;
  }
  head_b(outv, inv, vel, outv_b, inv_b, vel_b, N);
  return 0;
}
#endif*/

void head_b (double *u, double *u_1, double C, double D, double *u_b, double *u_1_b, int n){
  burgers1d_perf_b(u, u_b, u_1, u_1_b, D, C, n);
}
