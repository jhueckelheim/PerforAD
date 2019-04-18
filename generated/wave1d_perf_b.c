
#ifndef TAPENADE
#include <math.h>
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)

#define u(x) u[x]
#define u_b(x) u_b[x]
#define c(x) c[x]
#define u_1(x) u_1[x]
#define u_1_b(x) u_1_b[x]
#define u_2(x) u_2[x]
#define u_2_b(x) u_2_b[x]
void wave1d_perf_b(double* u, double* u_b, double* c, double* u_1, double* u_1_b, double* u_2, double* u_2_b, double D, int n) {
    int i;
    i=0;
    u_1_b(i) += D*c(i + 1)*u_b(i + 1);
    i=n - 2;
    u_1_b(i) += (-2*D*c(i) + 2.0)*u_b(i);
    u_2_b(i) += -u_b(i);
    u_1_b(i) += D*c(i - 1)*u_b(i - 1);
    i=1;
    u_1_b(i) += D*c(i + 1)*u_b(i + 1);
    u_1_b(i) += (-2*D*c(i) + 2.0)*u_b(i);
    u_2_b(i) += -u_b(i);
    i=n - 1;
    u_1_b(i) += D*c(i - 1)*u_b(i - 1);
    #pragma omp parallel for private(i)
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b(i) += D*c(i + 1)*u_b(i + 1);
        u_1_b(i) += (-2*D*c(i) + 2.0)*u_b(i);
        u_2_b(i) += -u_b(i);
        u_1_b(i) += D*c(i - 1)*u_b(i - 1);
    }
}