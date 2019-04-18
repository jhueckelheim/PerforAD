
#ifndef TAPENADE
#include <math.h>
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)

#define u(x,xx,xxx) u[x][xx][xxx]
#define c(x,xx,xxx) c[x][xx][xxx]
#define u_1(x,xx,xxx) u_1[x][xx][xxx]
#define u_2(x,xx,xxx) u_2[x][xx][xxx]
void wave3d(double*** u, double*** c, double*** u_1, double*** u_2, double D, int n) {
    int i;
    int j;
    int k;
    #pragma omp parallel for private(k,j,i)
    for ( i=1; i<=n - 2; i++ ) {
        for ( j=1; j<=n - 2; j++ ) {
            for ( k=1; k<=n - 2; k++ ) {
                u(i,j,k) += D*(-6*u_1(i, j, k) + u_1(i, j, k - 1) + u_1(i, j, k + 1) + u_1(i, j - 1, k) + u_1(i, j + 1, k) + u_1(i - 1, j, k) + u_1(i + 1, j, k))*c(i, j, k) + 2.0*u_1(i, j, k) - u_2(i, j, k);
            }
        }
    }
}