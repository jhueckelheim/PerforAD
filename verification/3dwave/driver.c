#include <stdio.h>
#include <stdlib.h>

#ifdef _CIVL
$input int N;
$input double	global_u[N][N][N];
$input double	global_c[N][N][N];
$input double global_u_1[N][N][N];
$input double global_u_2[N][N][N];
$input double global_D;
$input double 	global_u_b[N][N][N];
$output	double 	global_u_1_b[N][N][N];
$output double  global_u_2_b[N][N][N];
$assume(6<=N && N<=10);
#endif


void head_b(double*** u, double*** u_b, double*** c, double*** u_1, double*** u_1_b, double*** u_2, double*** u_2_b, double D, int n); 
double ***allocate3Ddouble(int d1, int d2, int d3){
 int i,j;
 double ***p = (double***) malloc(d1*sizeof(double**));
  for (i = 0; i < d1; i++){
    p[i] = (double**) malloc(d2*sizeof(double*));
    for (j = 0; j < d2; j++){
      p[i][j] = (double*) malloc(d3*sizeof(double));
    }
  }
  return p;
}

void free3Ddouble(double ***p, int d1, int d2, int d3){
 int i,j;
  for (i = 0; i < d1; i++){
    for (j = 0; j < d2; j++){
     free(p[i][j]);
    }
    free(p[i]);
  }
  free(p[i]);
}
/*
double ***allocate3Ddouble(int d1, int d2, int d3)
{
  int i,j, idx;
  double ***p = (double***) malloc(d1*sizeof(double**) + d1*d2*sizeof(double*) + d1 * d2 * d3 * sizeof(double));
  for (i = 0; i < d1; i++){
    p[i] = p+d1+d2*i;
  }
  for (i = 0; i < d1; i++){
    for (j = 0; j < d2; j++){
      p[i][j] = p + d1+d1*d2+d3*(i*d2+j);
    }
  }
  return p;
}

void free3Ddouble(double ***arr){
  free(arr);
}*/

int main()
{
  double  ***u, ***u_1, ***u_2, ***c;
  double D;
  double  ***u_b, ***u_1_b, ***u_2_b;
  int i,j,k;
  u = allocate3Ddouble(N,N,N);
  c = allocate3Ddouble(N,N,N);
  u_1 = allocate3Ddouble(N,N,N);
  u_2 = allocate3Ddouble(N,N,N);
  u_b = allocate3Ddouble(N,N,N);
  u_1_b = allocate3Ddouble(N,N,N);
  u_2_b = allocate3Ddouble(N,N,N);

  D = global_D;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        u[i][j][k] = global_u[i][j][k];
        c[i][j][k] = global_c[i][j][k];
        u_1[i][j][k] = global_u_1[i][j][k];
        u_2[i][j][k] = global_u_2[i][j][k];
        u_b[i][j][k] = global_u_b[i][j][k];
      }
    }
  }
  head_b(u, u_b, c, u_1, u_1_b, u_2, u_2_b, D, N);
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        global_u_1_b[i][j][k] = u_1_b[i][j][k];
        global_u_2_b[i][j][k] = u_2_b[i][j][k];
      }
    }
  }
  free3Ddouble(u,N,N,N);
  free3Ddouble(c,N,N,N);
  free3Ddouble(u_1,N,N,N);
  free3Ddouble(u_2,N,N,N);
  free3Ddouble(u_b,N,N,N);
  free3Ddouble(u_1_b,N,N,N);
  free3Ddouble(u_2_b,N,N,N);
  return 0;
}
