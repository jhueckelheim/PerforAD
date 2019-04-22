#ifndef TAPENADE
#include <math.h>
#else
double fmax(double a, double b) {
  if(a>b) return a;
  else return b;
}
double fmin(double a, double b) {
  if(a<b) return a;
  else return b;
}
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)

#define u(x) u[x]
#define u_1(x) u_1[x]
void burgers1d(double* u, double* u_1, double D, double C, int n) {
    int i;
    #pragma omp parallel for private(i)
    for ( i=1; i<=n - 2; i++ ) {
        u(i) += -C*((-u_1(i) + u_1(i + 1))*Min(0, u_1(i)) + (u_1(i) - u_1(i - 1))*Max(0, u_1(i))) + D*(-2.0*u_1(i) + u_1(i - 1) + u_1(i + 1)) + u_1(i);
    }
}
