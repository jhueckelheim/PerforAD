#include <stdio.h>

#include "ad_types.h"
#include "math.h"
#include "stdlib.h"

void head (double **outv, double **inv, double  **vel, int n);
void ad_head(DERIV_TYPE **outv,DERIV_TYPE **inv,DERIV_TYPE **vel,int n);

#ifdef DEBUG
void an_head (double **outv, double **inv,  double **vel, double **outv_b, double **inv_b, double **vel_b,int n);
#endif

double **allocate2Ddouble(int d1, int d2)
{
    int i;
    double *p = (double*) malloc(d1 * d2 * sizeof(double));
    double **q = (double**) malloc(d1 * sizeof(double*));
    for (i = 0; i < d1; i++)
    {
      int idx = d1 *i;
      q[i] = &p[idx];
    }
    return q;
} 

DERIV_TYPE **allocate2DDERIV_TYPE(int d1, int d2)
{
    int i;
    DERIV_TYPE *p = (DERIV_TYPE*) malloc(d1 * d2 * sizeof(DERIV_TYPE));
    DERIV_TYPE **q = (DERIV_TYPE**) malloc(d1 * sizeof(DERIV_TYPE*));
    for (i = 0; i < d1; i++)
    {
      int idx = d1 *i;
      q[i] = &p[idx];
    }
    return q;
}

void free2Ddouble(double **arr){
  free((arr[0]));
  free(arr);
}

void free2DDERIV_TYPE(DERIV_TYPE **arr){
  free((arr[0]));
  free(arr);
}


int main()
{
  double **outv, **inv, **vel;
  DERIV_TYPE **ad_outv, **ad_inv, **ad_vel;
  int i,j,k,l;
  outv = allocate2Ddouble(N,N);
  inv = allocate2Ddouble(N,N);
  vel = allocate2Ddouble(N,N);
  ad_outv = allocate2DDERIV_TYPE(N,N);
  ad_inv = allocate2DDERIV_TYPE(N,N);
  ad_vel = allocate2DDERIV_TYPE(N,N);
#ifdef DEBUG
  double **invph, **outvph;
  double h, **dd_inv;
  double **outv_b, **inv_b, **vel_b;
  invph = allocate2Ddouble(N,N);
  outvph = allocate2Ddouble(N,N);
  dd_inv = allocate2Ddouble(N,N);
  outv_b = allocate2Ddouble(N,N);
  inv_b = allocate2Ddouble(N,N);
  vel_b = allocate2Ddouble(N,N);

  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      outv[i][j] = 0.0;
      inv[i][j] = (i+1)*(j+1) * 0.5;
      invph[i][j] = (i+1)*(j+1) * 0.5;
      vel[i][j] = (i+1)*(j+1) * 0.5;
    }
  }

  /* dd step size */
  h=.0001;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
        for (l = 0; l < N; l++){
          outv[k][l] = 0.0;
          outvph[k][l] = 0.0;
        }
      } 
      invph[i][j] = inv[i][j] + h; 

      head(outv,inv, vel,N);
      head(outvph,invph, vel,N);
      dd_inv[i][j] = 0.0;
      for (k = 0; k < N; k++){
        for (l = 0; l < N; l++){
          dd_inv[i][j] += (outvph[k][l] - outv[k][l])/h;
        }
      }
      invph[i][j] = inv[i][j];
    }
  }
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      outv[i][j] = 0.0;
      inv[i][j] = (i+1)*(j+1) * 0.5;
      vel[i][j] = (i+1)*(j+1) * 0.5;
    }
  }
   
  head(outv,inv, vel,N);
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      outv_b[i][j] = 1.0;
      inv_b[i][j] = 0.0;
      vel_b[i][j] = 0.0;
    }
  }
  an_head(outv,inv, vel,outv_b,inv_b,vel_b, N);
#endif
  /* ----------------------------------------------------------
   * Compute dy/dx using AD 
   */
  ADIC_Init();
  __ADIC_TapeInit();
  
  
  // Set indpendent variables 
  ADIC_SetReverseMode();
  ADIC_SetDepMatrix(ad_outv, N, N);
  ADIC_SetIndepMatrix(ad_inv, N, N);
  ADIC_SetIndepMatrix(ad_vel, N, N);
  ADIC_SetIndepDone();
  // Initialize the value of the independent variable ad_x
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      inv[i][j] = (i+1)*(j+1) * 0.5;
      vel[i][j] = (i+1)*(j+1) * 0.5;
      DERIV_val(ad_outv[i][j]) =0.0;
      DERIV_val(ad_inv[i][j]) =inv[i][j];
      DERIV_val(ad_vel[i][j]) =vel[i][j];
    }
  }
  
  
  // Invoke AD function 
  our_rev_mode.tape = 1; 
  our_rev_mode.adjoint = 0; 
  ad_head(ad_outv,ad_inv, ad_vel,N);
#ifdef DEBUG
  printf("Primal Output\n");
  printf("i\tOrig\tADIC\n");
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      printf("%d %d\t%lf\t%lf \n", i,  j, outv[i][j], DERIV_val(ad_outv[i][j]));
    } 
  }
#endif
  our_rev_mode.tape = 0; 
  our_rev_mode.adjoint = 1; 
  ad_head(ad_outv,ad_inv,ad_vel,N);

#ifdef DEBUG
  printf("Adjoint Output inv \n");
  printf("i\tADIC\tPerforAD\tFD\n");
  double temp_adj;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      temp_adj = 0.0; 
      for (k = 0; k <ADIC_GRADVEC_LENGTH; k++) {
        temp_adj +=  DERIV_grad(ad_inv[i][j])[k];
      }
      printf("%d %d %lf\t%lf\t%lf \n", i, j,temp_adj, inv_b[i][j], dd_inv[i][j]); 
   }
  } 
  printf("Adjoint Output vel \n");
  printf("i\tADIC\tPerforAD\n");
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      temp_adj = 0.0; 
      for (k = 0; k <ADIC_GRADVEC_LENGTH; k++) {
        temp_adj +=  DERIV_grad(ad_vel[i][j])[k];
      }
      printf("%d %d %lf\t%lf \n", i, j, temp_adj, vel_b[i][j]); 
   }
  }
  free2Ddouble(invph);
  free2Ddouble(outvph);
  free2Ddouble(dd_inv);
  free2Ddouble(outv_b);
  free2Ddouble(inv_b);
  free2Ddouble(vel_b);
#endif  
  ADIC_Finalize();
  free2DDERIV_TYPE(ad_inv);
  free2DDERIV_TYPE(ad_vel);
  free2DDERIV_TYPE(ad_outv);
  free2Ddouble(inv);
  free2Ddouble(vel);
  free2Ddouble(outv);
  return 0;
}
#ifdef DEBUG
void an_head (double **outv, double **inv,  double **vel, double **outv_b, double **inv_b, double **vel_b,int n){
int i,j;
for ( j=1; j<=n - 2; j++ ) {
    i=0;
    inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
}
j=0;
i=n - 2;
inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
j=n - 2;
i=n - 2;
inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
j=1;
i=n - 2;
inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
j=n - 1;
i=n - 2;
inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
for ( j=2; j<=n - 3; j++ ) {
    i=n - 2;
    inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
    inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
    inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
    inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
}
j=0;
i=1;
inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
j=n - 2;
i=1;
inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
j=1;
i=1;
inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
j=n - 1;
i=1;
inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
for ( j=2; j<=n - 3; j++ ) {
    i=1;
    inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
    inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
    inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
    inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
}
for ( j=1; j<=n - 2; j++ ) {
    i=n - 1;
    inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
}
j=0;
for ( i=2; i<=n - 3; i++ ) {
    inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
}
j=n - 2;
for ( i=2; i<=n - 3; i++ ) {
    inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
    inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
    inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
    inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
}
j=1;
for ( i=2; i<=n - 3; i++ ) {
    inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
    inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
    inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
    inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
}
j=n - 1;
for ( i=2; i<=n - 3; i++ ) {
    inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
}
for ( j=2; j<=n - 3; j++ ) {
    for ( i=2; i<=n - 3; i++ ) {
        inv_b[i][j] += outv_b[i][ j + 1]*vel[i][ j + 1];
        inv_b[i][j] += outv_b[i + 1][ j]*vel[i + 1][ j];
        inv_b[i][j] += -4.0*outv_b[i][ j]*vel[i][ j];
        inv_b[i][j] += outv_b[i - 1][ j]*vel[i - 1][ j];
        inv_b[i][j] += outv_b[i][ j - 1]*vel[i][ j - 1];
    }
}
}
#endif

