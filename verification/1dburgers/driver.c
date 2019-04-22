#include <stdio.h>

#ifdef _CIVL
$input int N;
$input double	global_in_u[N];
$input	double 	global_in_u_1[N];
$input double   global_in_C;
$input  double  global_in_D;
$input	double 	global_u_b[N];
$output	double 	global_u_1_b[N];
$assume(6<=N && N<=10);
#endif

void head_b (double *u, double *u_1, double C, double D, double *u_b, double *u_1_b, int n);

int main()
{
  double u[N], u_1[N], C, D, u_b[N], u_1_b[N];
  int i;
  C = global_in_C;
  D = global_in_D;
  for (i = 0; i < N; i++){
    u[i] = global_in_u[i];
    u_1[i] = global_in_u_1[i];
    u_b[i] = global_u_b[i];
  }
  head_b(u, u_1, C, D, u_b, u_1_b, N);
  for (i = 0; i < N; i++){
    global_u_1_b[i] = u_1_b[i];
  }
  return 0;
}
