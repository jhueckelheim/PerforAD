#include "stdio.h"
void wave1d(double* u, double* u_1, double* u_2, double* c, double D, int n);
void wave1d_b(double *u, double *ub, double *u_1, double *u_1b, double *u_2, double *u_2b, double *c, double D, int n);
void wave1d_perf_b(double* u, double* u_1, double* u_2, double* c, double* u_b, double* u_1_b, double* u_2_b, double D, int n);

int main() {
  double u[N];
  double u_1[N];
  double u_2[N];
  double u_b[N];
  double u_1_b_Tapenade[N];
  double u_2_b_Tapenade[N];
  double u_1_b_PerforAD[N];
  double u_2_b_PerforAD[N];
  double c[N];
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
  wave1d(&u[0], &u_1[0], &u_2[0], &c[0], D, N);
  for(i=0; i<N; i++) {
    printf("%f ",u[i]);
  }
  printf("\n");

  // fill u_b with arbitrary values, initialise u_1_b and u_2_b with zeroes.
  for(i=0; i<N; i++) {
    u_1_b_Tapenade[i] = 0.0;
    u_2_b_Tapenade[i] = 0.0;
    u_1_b_PerforAD[i] = 0.0;
    u_2_b_PerforAD[i] = 0.0;
    u_b[i] = 2.0*0.001*i;
  }
  // call Tapenade adjoint
  wave1d_b(&u[0], &u_b[0], &u_1[0], &u_1_b_Tapenade[0], &u_2[0], &u_2_b_Tapenade[0], &c[0], D, N);
  // call PerforAD adjoint
  wave1d_perf_b(&u[0], &u_1[0], &u_2[0], &c[0], &u_b[0], &u_1_b_PerforAD[0], &u_2_b_PerforAD[0], D, N);
  // print results side by side
  for(i=0; i<N; i++) {
    printf("[%f, %f] [%f, %f]\n",u_1_b_Tapenade[i],u_1_b_PerforAD[i],u_2_b_Tapenade[i],u_2_b_PerforAD[i]);
  }
  printf("\n");
  return 0;
}
