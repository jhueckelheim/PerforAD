
#ifndef TAPENADE
#include <math.h>
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)

#define u(x) u[x]
#define u_1(x) u_1[x]
#define u_b(x) u_b[x]
#define u_1_b(x) u_1_b[x]
void burgers1d_perf_b(double* u, double* u_1, double* u_b, double* u_1_b, double C, double D, int n) {
    int i;
    i=0;
    u_1_b(i) += (C*Max(0, u_1(i + 1)) + D)*u_b(i + 1);
    i=n - 2;
    u_1_b(i) += (-C*((-u_1(i) + u_1(i + 1))*Heaviside(-u_1(i)) + (u_1(i) - u_1(i - 1))*Heaviside(u_1(i)) + Max(0, u_1(i)) - Min(0, u_1(i))) - 2.0*D + 1)*u_b(i);
    u_1_b(i) += (-C*Min(0, u_1(i - 1)) + D)*u_b(i - 1);
    i=1;
    u_1_b(i) += (C*Max(0, u_1(i + 1)) + D)*u_b(i + 1);
    u_1_b(i) += (-C*((-u_1(i) + u_1(i + 1))*Heaviside(-u_1(i)) + (u_1(i) - u_1(i - 1))*Heaviside(u_1(i)) + Max(0, u_1(i)) - Min(0, u_1(i))) - 2.0*D + 1)*u_b(i);
    i=n - 1;
    u_1_b(i) += (-C*Min(0, u_1(i - 1)) + D)*u_b(i - 1);
    #pragma omp for private(i)
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b(i) += (C*Max(0, u_1(i + 1)) + D)*u_b(i + 1);
        u_1_b(i) += (-C*((-u_1(i) + u_1(i + 1))*Heaviside(-u_1(i)) + (u_1(i) - u_1(i - 1))*Heaviside(u_1(i)) + Max(0, u_1(i)) - Min(0, u_1(i))) - 2.0*D + 1)*u_b(i);
        u_1_b(i) += (-C*Min(0, u_1(i - 1)) + D)*u_b(i - 1);
    }
}