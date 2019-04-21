#include "stdio.h"
#include <omp.h>


void wave1d(double* u, double* u_1, double* u_2, double* c, double D, int n);
void wave1d_b(double *u, double *ub, double *u_1, double *u_1b, double *u_2, double *u_2b, double *c, double D, int n);
void wave1d_perf_b(double* u, double* u_b, double* c, double* u_1, double* u_1_b, double* u_2, double* u_2_b, double D, int n);

int main(int argc, char *argv[]) {
  int N;
  if(argc==0) {
    N = 500000000;
  } else {
    N = atoi(argv[1]);
  }
  double *u = malloc(N*sizeof(double));
  double *u_1 = malloc(N*sizeof(double));
  double *u_2 = malloc(N*sizeof(double));
  double *u_b = malloc(N*sizeof(double));
  double *u_1_b_Tapenade = malloc(N*sizeof(double));
  double *u_2_b_Tapenade = malloc(N*sizeof(double));
  double *u_1_b_PerforAD = malloc(N*sizeof(double));
  double *u_2_b_PerforAD = malloc(N*sizeof(double));
  double *c = malloc(N*sizeof(double));
  int i;
  double D;

  // fill inputs with arbitrary values.
  D = 3.1415;
  for(i=0; i<N; i++) {
    u_1[i] = 0.2*i*i;
    u_2[i] = -0.7*i*i;
    c[i] = 2.0*0.001*i;
  }
  // run primal and print its output.
  double start_time = omp_get_wtime();
  wave1d(&u[0], &u_1[0], &u_2[0], &c[0], D, N);
  double time = omp_get_wtime() - start_time;
  printf("Forward: %f\n", time);
  //for(i=0; i<N; i++) {
  //  printf("%f ",u[i]);
  //}
  //printf("\n");

  // fill u_b with arbitrary values, initialise u_1_b and u_2_b with zeroes.
  for(i=0; i<N; i++) {
    u_1_b_Tapenade[i] = 0.0;
    u_2_b_Tapenade[i] = 0.0;
    u_1_b_PerforAD[i] = 0.0;
    u_2_b_PerforAD[i] = 0.0;
    u_b[i] = 2.0*0.001*i;
  }
  // call Tapenade adjoint
  start_time = omp_get_wtime();
  wave1d_b(&u[0], &u_b[0], &u_1[0], &u_1_b_Tapenade[0], &u_2[0], &u_2_b_Tapenade[0], &c[0], D, N);
  time = omp_get_wtime() - start_time;
  printf("Tapenade: %f\n", time);
  // call PerforAD adjoint
  start_time = omp_get_wtime();
  wave1d_perf_b(&u[0], &u_b[0], &c[0], &u_1[0], &u_1_b_PerforAD[0], &u_2[0], &u_2_b_PerforAD[0], D, N);
  time = omp_get_wtime() - start_time;
  printf("PerforAD: %f\n", time);
  // print results side by side
  //for(i=0; i<N; i++) {
  //  printf("[%f, %f] [%f, %f]\n",u_1_b_Tapenade[i],u_1_b_PerforAD[i],u_2_b_Tapenade[i],u_2_b_PerforAD[i]);
  //}
  //printf("\n");
  free(u);
  free(u_b);
  free(c);
  free(u_1);
  free(u_1_b_Tapenade);
  free(u_2);
  free(u_2_b_Tapenade);
  free(u_1_b_PerforAD);
  free(u_2_b_PerforAD);

  return 0;
}
