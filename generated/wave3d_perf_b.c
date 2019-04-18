
#ifndef TAPENADE
#include <math.h>
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)

#define u(x,xx,xxx) u[x][xx][xxx]
#define u_1(x,xx,xxx) u_1[x][xx][xxx]
#define u_2(x,xx,xxx) u_2[x][xx][xxx]
#define c(x,xx,xxx) c[x][xx][xxx]
#define u_b(x,xx,xxx) u_b[x][xx][xxx]
#define u_1_b(x,xx,xxx) u_1_b[x][xx][xxx]
#define u_2_b(x,xx,xxx) u_2_b[x][xx][xxx]
void wave3d_perf_b(double*** u, double*** u_1, double*** u_2, double*** c, double*** u_b, double*** u_1_b, double*** u_2_b, double D, int n) {
    int i;
    int j;
    int k;
    i=0;
    #pragma omp parallel for private(k,j)
    for ( j=1; j<=n - 2; j++ ) {
        for ( k=1; k<=n - 2; k++ ) {
            u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        }
    }
    i=n - 2;
    j=0;
    #pragma omp parallel for private(k)
    for ( k=1; k<=n - 2; k++ ) {
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
    }
    i=n - 2;
    j=n - 2;
    k=0;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    i=n - 2;
    j=n - 2;
    k=n - 2;
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
    u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=n - 2;
    j=n - 2;
    k=1;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
    u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    i=n - 2;
    j=n - 2;
    k=n - 1;
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=n - 2;
    j=n - 2;
    #pragma omp parallel for private(k)
    for ( k=2; k<=n - 3; k++ ) {
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=n - 2;
    j=1;
    k=0;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    i=n - 2;
    j=1;
    k=n - 2;
    u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=n - 2;
    j=1;
    k=1;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
    i=n - 2;
    j=1;
    k=n - 1;
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=n - 2;
    j=1;
    #pragma omp parallel for private(k)
    for ( k=2; k<=n - 3; k++ ) {
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=n - 2;
    j=n - 1;
    #pragma omp parallel for private(k)
    for ( k=1; k<=n - 2; k++ ) {
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    }
    i=n - 2;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=0;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    }
    i=n - 2;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=n - 2;
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=n - 2;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=1;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    }
    i=n - 2;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=n - 1;
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=n - 2;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        for ( k=2; k<=n - 3; k++ ) {
            u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
            u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
            u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
            u_2_b(i,j,k) += -u_b(i, j, k);
            u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
            u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
            u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
        }
    }
    i=1;
    j=0;
    #pragma omp parallel for private(k)
    for ( k=1; k<=n - 2; k++ ) {
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
    }
    i=1;
    j=n - 2;
    k=0;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    i=1;
    j=n - 2;
    k=n - 2;
    u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=1;
    j=n - 2;
    k=1;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    i=1;
    j=n - 2;
    k=n - 1;
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=1;
    j=n - 2;
    #pragma omp parallel for private(k)
    for ( k=2; k<=n - 3; k++ ) {
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=1;
    j=1;
    k=0;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    i=1;
    j=1;
    k=n - 2;
    u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
    u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=1;
    j=1;
    k=1;
    u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
    u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
    u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
    u_2_b(i,j,k) += -u_b(i, j, k);
    i=1;
    j=1;
    k=n - 1;
    u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    i=1;
    j=1;
    #pragma omp parallel for private(k)
    for ( k=2; k<=n - 3; k++ ) {
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=1;
    j=n - 1;
    #pragma omp parallel for private(k)
    for ( k=1; k<=n - 2; k++ ) {
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    }
    i=1;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=0;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    }
    i=1;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=n - 2;
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=1;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=1;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    }
    i=1;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        k=n - 1;
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    i=1;
    #pragma omp parallel for private(k,j)
    for ( j=2; j<=n - 3; j++ ) {
        for ( k=2; k<=n - 3; k++ ) {
            u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
            u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
            u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
            u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
            u_2_b(i,j,k) += -u_b(i, j, k);
            u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
            u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
        }
    }
    i=n - 1;
    #pragma omp parallel for private(k,j)
    for ( j=1; j<=n - 2; j++ ) {
        for ( k=1; k<=n - 2; k++ ) {
            u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=0;
        for ( k=1; k<=n - 2; k++ ) {
            u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=n - 2;
        k=0;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=n - 2;
        k=n - 2;
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=n - 2;
        k=1;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=n - 2;
        k=n - 1;
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=n - 2;
        for ( k=2; k<=n - 3; k++ ) {
            u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
            u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
            u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
            u_2_b(i,j,k) += -u_b(i, j, k);
            u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
            u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
            u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=1;
        k=0;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=1;
        k=n - 2;
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=1;
        k=1;
        u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
        u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
        u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
        u_2_b(i,j,k) += -u_b(i, j, k);
        u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=1;
        k=n - 1;
        u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=1;
        for ( k=2; k<=n - 3; k++ ) {
            u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
            u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
            u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
            u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
            u_2_b(i,j,k) += -u_b(i, j, k);
            u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
            u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        j=n - 1;
        for ( k=1; k<=n - 2; k++ ) {
            u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        for ( j=2; j<=n - 3; j++ ) {
            k=0;
            u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        for ( j=2; j<=n - 3; j++ ) {
            k=n - 2;
            u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
            u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
            u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
            u_2_b(i,j,k) += -u_b(i, j, k);
            u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
            u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
            u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        for ( j=2; j<=n - 3; j++ ) {
            k=1;
            u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
            u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
            u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
            u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
            u_2_b(i,j,k) += -u_b(i, j, k);
            u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
            u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        for ( j=2; j<=n - 3; j++ ) {
            k=n - 1;
            u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
        }
    }
    #pragma omp parallel for private(k,j,i)
    for ( i=2; i<=n - 3; i++ ) {
        for ( j=2; j<=n - 3; j++ ) {
            for ( k=2; k<=n - 3; k++ ) {
                u_1_b(i,j,k) += D*c(i, j, k + 1)*u_b(i, j, k + 1);
                u_1_b(i,j,k) += D*c(i, j + 1, k)*u_b(i, j + 1, k);
                u_1_b(i,j,k) += D*c(i + 1, j, k)*u_b(i + 1, j, k);
                u_1_b(i,j,k) += (-6*D*c(i, j, k) + 2.0)*u_b(i, j, k);
                u_2_b(i,j,k) += -u_b(i, j, k);
                u_1_b(i,j,k) += D*c(i - 1, j, k)*u_b(i - 1, j, k);
                u_1_b(i,j,k) += D*c(i, j - 1, k)*u_b(i, j - 1, k);
                u_1_b(i,j,k) += D*c(i, j, k - 1)*u_b(i, j, k - 1);
            }
        }
    }
}