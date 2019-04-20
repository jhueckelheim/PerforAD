#include <time.h>
#include <stdlib.h>
#include <omp.h>


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
  
  init_random(u, n);
  init_random(u_1, n);

  double C = 1.2;
  double D = 3.23;
  
  
  double start_time = omp_get_wtime();
  burgers1d(u, u_1, C, D, n);
  double time = omp_get_wtime() - start_time;
  printf("%f", time);
}
