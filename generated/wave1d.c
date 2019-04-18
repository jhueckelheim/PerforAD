
#ifndef TAPENADE
#include <math.h>
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)

#define u(x) u[x]
#define u_1(x) u_1[x]
#define u_2(x) u_2[x]
#define c(x) c[x]
void wave1d(double* u, double* u_1, double* u_2, double* c, double D, int n) {
    int i;
    #pragma omp parallel for private(i)
    for ( i=1; i<=n - 2; i++ ) {
        u(i) += D*(-2*u_1(i) + u_1(i - 1) + u_1(i + 1))*c(i) + 2.0*u_1(i) - u_2(i);
    }
}