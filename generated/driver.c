#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void burgers1d(double* u, double* u_1, double D, double C, int n);

void burgers1d_b(double *u, double *ub, double *u_1, double *u_1b, double D,double C, int n);

void burgers1d_perf_b(double* u, double* u_b, double* u_1, double* u_1_b, double D, double C, int n);

void init_random(double *array, int n) {
  srand(time(NULL));
  int i;
  for(i=0; i<n; i++) {
    array[i] = rand();
  }
}

void main(int argc, char *argv[]) {
  int n;
  if(argc==0) {
    n = 500000000;
  } else {
    n = atoi(argv[1]);
  }
  double *u;
  u = malloc(n*sizeof(double));
  double *ub =  malloc(n*sizeof(double));
  double *u_1 = malloc(n*sizeof(double));
  double *u_1b = malloc(n*sizeof(double));

  init_random(u, n);
  init_random(u_1, n);

  init_random(ub, n);
  init_random(u_1b, n);

  double C = 1.2;
  double D = 3.23;
  
  
  double start_time = omp_get_wtime();
  burgers1d(u, u_1, C, D, n);
  double time = omp_get_wtime() - start_time;
  printf("Forward: %f\n", time);

  start_time = omp_get_wtime();
  burgers1d_b(u, ub, u_1, u_1b, D, C, n);
  time = omp_get_wtime() - start_time;

  printf("Tapenade: %f\n", time);

  start_time = omp_get_wtime();
  burgers1d_perf_b(u, ub, u_1, u_1b, D, C, n);
  time = omp_get_wtime() - start_time;

  printf("PerforAd: %f\n", time);
  free(u);
  free(u_1);
  free(ub);
  free(u_1b);
}
