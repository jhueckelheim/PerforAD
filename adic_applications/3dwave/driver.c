#include <stdio.h>

#include "ad_types.h"
#include "math.h"
#include "stdlib.h"

void head (double ***u, double ***u_1, double  ***u_2, double ***c, double D, int n);
void ad_head(DERIV_TYPE ***u,DERIV_TYPE ***u_1,DERIV_TYPE ***u_2,DERIV_TYPE ***c,DERIV_TYPE *D,int n);

#ifdef DEBUG
void an_head (double ***u, double ***u_1, double  ***u_2, double ***c, double ***u_b, double ***u_1_b, double  ***u_2_b, double ***c_b, double D, int n);
#endif

double ***allocate3Ddouble(int d1, int d2, int d3)
{
    int i,j, idx;
    double ***p = (double***) malloc(d1*sizeof(double**) + d1*d2*sizeof(double*) + d1 * d2 * d3 * sizeof(double));
    //printf("D &p:%ld p:%ld sizeof(double): %ld \n", &p, p, sizeof(double)); 
    for (i = 0; i < d1; i++){
      p[i] = p+d1+d2*i;
      //printf("D &p[%d]:%ld p[%d]:%ld p+d1+i:%ld \n", i, &(p[i]),i, p[i], p+d1+d2*i); 
    }
    for (i = 0; i < d1; i++){
      for (j = 0; j < d2; j++){
        p[i][j] = p + d1+d1*d2+d3*(i*d2+j);
        //printf("D &(p[%d][%d]): %ld p[%d][%d]: %ld  d3*sizeof(double)*(i*d2+j): %ld \n",i,j, &(p[i][j]), i,j, p[i][j], p+ d1+d1*d2+d3*sizeof(double)*(i*d2+j));
      }
    }
    return p;
} 

DERIV_TYPE ***allocate3DDERIV_TYPE(int d1, int d2, int d3)
{
    int i,j, idx;
    DERIV_TYPE ***p = (DERIV_TYPE***) malloc(d1*sizeof(DERIV_TYPE**) + d1*d2*sizeof(DERIV_TYPE*) + d1 * d2 * d3 * sizeof(DERIV_TYPE));
    //printf("D &p:%ld p:%ld sizeof(DERIV_TYPE): %ld \n", &p, p, sizeof(DERIV_TYPE)); 
    for (i = 0; i < d1; i++){
      p[i] = p+d1+d2*i;
      //printf("D &p[%d]:%ld p[%d]:%ld p+d1+i:%ld \n", i, &(p[i]),i, p[i], p+d1+d2*i); 
    }
    for (i = 0; i < d1; i++){
      for (j = 0; j < d2; j++){
        p[i][j] = p + d1+d1*d2+d3*sizeof(DERIV_TYPE)/8*(i*d2+j);
        //printf("D &(p[%d][%d]): %ld p[%d][%d]: %ld  d3*sizeof(DERIV_TYPE)*(i*d2+j): %ld \n",i,j, &(p[i][j]), i,j, p[i][j], p+ d1+d1*d2+d3*sizeof(DERIV_TYPE)*(i*d2+j));
      }
    }
    return p;
}
/* 
DERIV_TYPE ***allocate3Ddouble(int d1, int d2, int d3)
{
    int i,j, idx;
    double ***q = (double***) malloc(d1 * sizeof(double**));
    for (i = 0; i < d1; i++){
      q[i] = (double**) malloc(d2 * sizeof(double*));
      for (j = 0; j < d2; j++){
        q[i][j] = (double*) malloc(d3 * sizeof(double));
      }
    }
    return q;
} 


DERIV_TYPE ***allocate3DDERIV_TYPE(int d1, int d2, int d3)
{
    int i,j, idx;
    DERIV_TYPE ***q = (DERIV_TYPE***) malloc(d1 * sizeof(DERIV_TYPE**));
    for (i = 0; i < d1; i++){
      q[i] = (DERIV_TYPE**) malloc(d2 * sizeof(DERIV_TYPE*));
      for (j = 0; j < d2; j++){
        q[i][j] = (DERIV_TYPE*) malloc(d3 * sizeof(DERIV_TYPE));
      }
    }
    return q;
} */

void free3Ddouble(double ***arr){
  free(arr);
}

void free3DDERIV_TYPE(DERIV_TYPE ***arr){
  free(arr);
}

int main()
{
  double ***u, ***u_1, ***u_2, ***c, D;
  DERIV_TYPE ***ad_u, ***ad_u_1, ***ad_u_2, ***ad_c, ad_D;
  int i,j,k,l,m,n;
  u = allocate3Ddouble(N,N,N);
  u_1 = allocate3Ddouble(N,N,N);
  u_2 = allocate3Ddouble(N,N,N);
  c = allocate3Ddouble(N,N,N);
  ad_u = allocate3DDERIV_TYPE(N,N,N);
  ad_u_1 = allocate3DDERIV_TYPE(N,N,N);
  ad_u_2 = allocate3DDERIV_TYPE(N,N,N);
  ad_c = allocate3DDERIV_TYPE(N,N,N);
#ifdef DEBUG
  double ***u_1ph, ***uph;
  double h, ***dd_u_1;
  double ***u_b, ***u_1_b, ***u_2_b, ***c_b;
  u_1ph = allocate3Ddouble(N,N,N);
  uph = allocate3Ddouble(N,N,N);
  dd_u_1 = allocate3Ddouble(N,N,N);
  u_b = allocate3Ddouble(N,N,N);
  u_1_b = allocate3Ddouble(N,N,N);
  u_2_b = allocate3Ddouble(N,N,N);
  c_b = allocate3Ddouble(N,N,N);

  D = 0.5;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        u[i][j][k] = 0.0;
        u_1[i][j][k] = (i+1)*(j+1) * 0.5;
        u_1ph[i][j][k] = (i+1)*(j+1) * 0.5;
        u_2[i][j][k] = (i+1)*(j+1) * 0.5;
        c[i][j][k] = (i+1)*(j+1) * 0.5;
      }
    }
  }

  /* dd step size */
  h=.0000001;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        for (l = 0; l < N; l++){
          for (m = 0; m < N; m++){
            for (n = 0; n < N; n++){
              u[l][m][n] = 0.0;
              uph[l][m][n] = 0.0;
            }
          }
        } 
        u_1ph[i][j][k] = u_1[i][j][k] + h; 

        head(u, u_1, u_2, c, D, N);
        head(uph, u_1ph, u_2, c, D, N);
        dd_u_1[i][j][k] = 0.0;
        for (l = 0; l < N; l++){
          for (m = 0; m < N; m++){
            for (n = 0; n < N; n++){
              dd_u_1[i][j][k] += (uph[l][m][n] - u[l][m][n])/h;
            }
          }
        }
        u_1ph[i][j][k] = u_1[i][j][k];
      }
    }
  }
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        u[i][j][k] = 0.0;
        u_1[i][j][k] = (i+1)*(j+1) * 0.5;
        u_2[i][j][k] = (i+1)*(j+1) * 0.5;
        c[i][j][k] = (i+1)*(j+1) * 0.5;
      }
    }
  }

  head(u, u_1, u_2, c, D, N);
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        u_b[i][j][k] = 1.0;
        u_1_b[i][j][k] = 0.0;
        u_2_b[i][j][k] = 0.0;
        c_b[i][j][k] = 0.0;
      }
    }
  }
  an_head(u, u_1, u_2, c, u_b, u_1_b, u_2_b, c_b, D, N);
#endif
  /* ----------------------------------------------------------
   * Compute dy/dx using AD 
   */
  ADIC_Init();
  __ADIC_TapeInit();
  
  
  // Set indpendent variables 
  ADIC_SetReverseMode();
  ADIC_SetDep3D(ad_u, N, N,N);
  ADIC_SetIndep3D(ad_u_1, N, N, N);
  ADIC_SetIndep3D(ad_u_2, N, N, N);
  ADIC_SetIndep3D(ad_c, N, N, N);
  ADIC_SetIndepDone();
  // Initialize parameters
  DERIV_val(ad_D) = 0.5;
  ZeroDeriv(ad_D);
  ZeroDerivthree_tensor(ad_c, N, N, N);
  // Initialize the value of the independent variable ad_x
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        u_1[i][j][k] = (i+1)*(j+1) * 0.5;
        u_2[i][j][k] = (i+1)*(j+1) * 0.5;
        c[i][j][k] = (i+1)*(j+1) * 0.5;
        DERIV_val(ad_u[i][j][k]) =0.0;
        DERIV_val(ad_u_1[i][j][k]) =u_1[i][j][k];
        DERIV_val(ad_u_2[i][j][k]) =u_2[i][j][k];
        DERIV_val(ad_c[i][j][k]) =c[i][j][k];
      }
    }
  }
  
  
  // Invoke AD function 
  our_rev_mode.tape = 1; 
  our_rev_mode.adjoint = 0; 
  ad_head(ad_u,ad_u_1, ad_u_2,ad_c,&D,N);
#ifdef DEBUG
  printf("Primal Output\n");
  printf("i\tOrig\tADIC\n");
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        printf("%d %d %d\t%lf\t%lf \n", i,  j, k, u[i][j][k], DERIV_val(ad_u[i][j][k]));
      }
    } 
  }
#endif
  our_rev_mode.tape = 0; 
  our_rev_mode.adjoint = 1; 
  ad_head(ad_u,ad_u_1, ad_u_2,ad_c,&D,N);

#ifdef DEBUG
  printf("Adjoint Output u_1 \n");
  printf("i\tADIC\tPerforAD\tFD\n");
  double temp_adj;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++) {
        temp_adj = 0.0; 
        for (l = 0; l <ADIC_GRADVEC_LENGTH; l++) {
          temp_adj +=  DERIV_grad(ad_u_1[i][j][k])[l];
        }
        printf("%d %d %d %lf\t%lf\t%lf \n", i, j,k, temp_adj, u_1_b[i][j][k], dd_u_1[i][j][k]); 
      }
    }
  } 
  printf("Adjoint Output u_2 \n");
  printf("i\tADIC\tPerforAD\n");
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++) {
        temp_adj = 0.0; 
        for (l = 0; l <ADIC_GRADVEC_LENGTH; l++) {
          temp_adj +=  DERIV_grad(ad_u_2[i][j][k])[l];
        }
        printf("%d %d %d %lf\t%lf \n", i, j, k, temp_adj, u_2_b[i][j][k]); 
      }
    }
  } 
  free3Ddouble(u_1ph);
  free3Ddouble(uph);
  free3Ddouble(dd_u_1);
  free3Ddouble(u_b);
  free3Ddouble(u_1_b);
  free3Ddouble(u_2_b);
  free3Ddouble(c_b);
#endif  
  ADIC_Finalize();
  free3DDERIV_TYPE(ad_u_1);
  free3DDERIV_TYPE(ad_u_2);
  free3DDERIV_TYPE(ad_u);
  free3DDERIV_TYPE(ad_c);
  free3Ddouble(u_1);
  free3Ddouble(u_2);
  free3Ddouble(u);
  free3Ddouble(c);
  return 0;
}
#ifdef DEBUG
void an_head (double ***u, double ***u_1, double  ***u_2, double ***c, double ***u_b, double ***u_1_b, double  ***u_2_b, double ***c_b, double D, int n){
int i,j,k;
for ( k=1; k<=n - 2; k++ ) {
    for ( j=1; j<=n - 2; j++ ) {
        i=0;
        u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    }
}
for ( k=1; k<=n - 2; k++ ) {
    j=0;
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
}
k=0;
j=n - 2;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
k=n - 2;
j=n - 2;
i=n - 2;
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
k=1;
j=n - 2;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
k=n - 1;
j=n - 2;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
for ( k=2; k<=n - 3; k++ ) {
    j=n - 2;
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
k=0;
j=1;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
k=n - 2;
j=1;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
k=1;
j=1;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
k=n - 1;
j=1;
i=n - 2;
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
for ( k=2; k<=n - 3; k++ ) {
    j=1;
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
for ( k=1; k<=n - 2; k++ ) {
    j=n - 1;
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
}
k=0;
for ( j=2; j<=n - 3; j++ ) {
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
}
k=n - 2;
for ( j=2; j<=n - 3; j++ ) {
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
k=1;
for ( j=2; j<=n - 3; j++ ) {
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
}
k=n - 1;
for ( j=2; j<=n - 3; j++ ) {
    i=n - 2;
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
for ( k=2; k<=n - 3; k++ ) {
    for ( j=2; j<=n - 3; j++ ) {
        i=n - 2;
        u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
        u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
        u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
        u_2_b[i][j][k] += -u_b[i][ j][ k];
        u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
        u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
        u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
    }
}
for ( k=1; k<=n - 2; k++ ) {
    j=0;
    i=1;
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
}
k=0;
j=n - 2;
i=1;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
k=n - 2;
j=n - 2;
i=1;
u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
k=1;
j=n - 2;
i=1;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
k=n - 1;
j=n - 2;
i=1;
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
for ( k=2; k<=n - 3; k++ ) {
    j=n - 2;
    i=1;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
k=0;
j=1;
i=1;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
k=n - 2;
j=1;
i=1;
u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
k=1;
j=1;
i=1;
u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
u_2_b[i][j][k] += -u_b[i][ j][ k];
k=n - 1;
j=1;
i=1;
u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
for ( k=2; k<=n - 3; k++ ) {
    j=1;
    i=1;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
for ( k=1; k<=n - 2; k++ ) {
    j=n - 1;
    i=1;
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
}
k=0;
for ( j=2; j<=n - 3; j++ ) {
    i=1;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
}
k=n - 2;
for ( j=2; j<=n - 3; j++ ) {
    i=1;
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
k=1;
for ( j=2; j<=n - 3; j++ ) {
    i=1;
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
}
k=n - 1;
for ( j=2; j<=n - 3; j++ ) {
    i=1;
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
for ( k=2; k<=n - 3; k++ ) {
    for ( j=2; j<=n - 3; j++ ) {
        i=1;
        u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
        u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
        u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
        u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
        u_2_b[i][j][k] += -u_b[i][ j][ k];
        u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
        u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
    }
}
for ( k=1; k<=n - 2; k++ ) {
    for ( j=1; j<=n - 2; j++ ) {
        i=n - 1;
        u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    }
}
for ( k=1; k<=n - 2; k++ ) {
    j=0;
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    }
}
k=0;
j=n - 2;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
}
k=n - 2;
j=n - 2;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
k=1;
j=n - 2;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
}
k=n - 1;
j=n - 2;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
for ( k=2; k<=n - 3; k++ ) {
    j=n - 2;
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
        u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
        u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
        u_2_b[i][j][k] += -u_b[i][ j][ k];
        u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
        u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
        u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
    }
}
k=0;
j=1;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
}
k=n - 2;
j=1;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
k=1;
j=1;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
    u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
    u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
    u_2_b[i][j][k] += -u_b[i][ j][ k];
    u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
}
k=n - 1;
j=1;
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
}
for ( k=2; k<=n - 3; k++ ) {
    j=1;
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
        u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
        u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
        u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
        u_2_b[i][j][k] += -u_b[i][ j][ k];
        u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
        u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
    }
}
for ( k=1; k<=n - 2; k++ ) {
    j=n - 1;
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    }
}
k=0;
for ( j=2; j<=n - 3; j++ ) {
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
    }
}
k=n - 2;
for ( j=2; j<=n - 3; j++ ) {
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
        u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
        u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
        u_2_b[i][j][k] += -u_b[i][ j][ k];
        u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
        u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
        u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
    }
}
k=1;
for ( j=2; j<=n - 3; j++ ) {
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
        u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
        u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
        u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
        u_2_b[i][j][k] += -u_b[i][ j][ k];
        u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
        u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
    }
}
k=n - 1;
for ( j=2; j<=n - 3; j++ ) {
    for ( i=2; i<=n - 3; i++ ) {
        u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
    }
}
for ( k=2; k<=n - 3; k++ ) {
    for ( j=2; j<=n - 3; j++ ) {
        for ( i=2; i<=n - 3; i++ ) {
            u_1_b[i][j][k] += D*c[i][ j][ k + 1]*u_b[i][ j][ k + 1];
            u_1_b[i][j][k] += D*c[i][ j + 1][ k]*u_b[i][ j + 1][ k];
            u_1_b[i][j][k] += D*c[i + 1][ j][ k]*u_b[i + 1][ j][ k];
            u_1_b[i][j][k] += (-6*D*c[i][ j][ k] + 2.0)*u_b[i][ j][ k];
            u_2_b[i][j][k] += -u_b[i][ j][ k];
            u_1_b[i][j][k] += D*c[i - 1][ j][ k]*u_b[i - 1][ j][ k];
            u_1_b[i][j][k] += D*c[i][ j - 1][ k]*u_b[i][ j - 1][ k];
            u_1_b[i][j][k] += D*c[i][ j][ k - 1]*u_b[i][ j][ k - 1];
        }
    }
}
}
#endif
