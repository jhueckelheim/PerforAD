#include "stdio.h"
#include <omp.h>
#include <time.h>

void wave3d(double* u, double* c, double* u_1, double* u_2, double D, int n);
void wave3d_b(double *u, double *ub, double *u_1, double *u_1b, double *u_2, double *u_2b, double *c, double D, int n);
void wave3d_perf_b(double* u, double* u_b, double* c, double* u_1, double* u_1_b, double* u_2, double* u_2_b, double D, int n);

void init_random(double *vec, int n) {
  double (*array)[n][n] = (double (*)[n][n]) vec;
  srand(time(NULL));
  int i, j, k;
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      for(k=0; k<n; k++) {
	array[i][j][k] = rand();
      }
    }
  }
}

void main(int argc, char *argv[]) {
  int n;
  if(argc==0) {
    n = 100;
  } else {
    n = atoi(argv[1]);
  }
  double (*u) = malloc(sizeof(double)*n*n*n);
  double (*c) = malloc(sizeof(double)*n*n*n);
  double (*u_1) = malloc(sizeof(double)*n*n*n);
  double (*u_2) = malloc(sizeof(double)*n*n*n);
  double (*ub) = malloc(sizeof(double)*n*n*n);
  double (*u_1b) = malloc(sizeof(double)*n*n*n);
  double (*u_2b) = malloc(sizeof(double)*n*n*n);
 
  init_random(u, n);
  init_random(c, n);
  init_random(u_1, n);
  init_random(u_2, n);

  init_random(ub, n);
  init_random(u_1b, n);
  init_random(u_2b, n);

  double C = 1.2;
  double D = 3.23;
  
  
  double start_time = omp_get_wtime();
  wave3d(u, c, u_1, u_2, D, n);
  double time = omp_get_wtime() - start_time;
  printf("Forward: %f\n", time);

  start_time = omp_get_wtime();
  wave3d_b(u, ub, u_1, u_1b, u_2, u_2b, c, D, n);
  time = omp_get_wtime() - start_time;

  printf("Tapenade: %f\n", time);

  start_time = omp_get_wtime();
  wave3d_perf_b(u, ub, u_1, u_1b, u_2, u_2b, c, D, n);
  time = omp_get_wtime() - start_time;

  printf("PerforAd: %f\n", time);
  free(u);
  free(u_1);
  free(c);
  free(u_2);
  free(ub);
  free(u_1b);
  free(u_2b);
}
