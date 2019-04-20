#include <stdio.h>
#include "ad_types.h"

void ad_head(DERIV_TYPE *outv,DERIV_TYPE *invec,DERIV_TYPE *vel,int n);

void head_b (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b);
#ifndef _CIVL
int main()
{
  double outv[N], inv[N], vel[N];
  double outv_b[N], inv_b[N], vel_b[N];
  int i,j;
  for (i = 0; i < N; i++){
    outv[i] = 0.0;
    inv[i] = (i+1)*(i+1) * 0.5;
    vel[i] = (i+1)*(i+1) * 0.5;
    outv_b[i] = 1.0;
    inv_b[i] = 0.0;
    vel_b[i] = 0.0;
  }
  head_b(outv, inv, vel, outv_b, inv_b, vel_b);
  return 0;
}
#endif

void head_b (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b)
{
  DERIV_TYPE *ad_outv, *ad_inv, *ad_vel;
  int i,j;

  ad_outv = (DERIV_TYPE *)malloc(N*sizeof(DERIV_TYPE));
  ad_inv = (DERIV_TYPE *)malloc(N*sizeof(DERIV_TYPE));
  ad_vel = (DERIV_TYPE *)malloc(N*sizeof(DERIV_TYPE));
  /* ----------------------------------------------------------
   * Compute dy/dx using AD 
   */
  ADIC_Init();
  __ADIC_TapeInit();
  
  
  // Set indpendent variables 
  ADIC_SetReverseMode();
  ADIC_SetDepArray(ad_outv, N);
  ADIC_SetIndepArray(ad_inv, N);
  ADIC_SetIndepArray(ad_vel, N);
  ADIC_SetIndepDone();
  // Initialize the value of the independent variable ad_x
  for (i = 0; i < N; i++){
    DERIV_val(ad_outv[i]) = outv[i];
    DERIV_val(ad_inv[i]) =inv[i];
    DERIV_val(ad_vel[i]) =vel[i];
  }
  
  // Invoke AD function 
  our_rev_mode.tape = 1; 
  our_rev_mode.adjoint = 0; 
  ad_head(ad_outv,ad_inv, ad_vel,N);
#ifdef DEBUG
  printf("Primal Output\n");
  printf("i\tOrig\tADIC\n");
  for (i = 0; i < N; i++){
    printf("%d\t%lf \n", i, DERIV_val(ad_outv[i]));
   } 
#endif
  our_rev_mode.tape = 0; 
  our_rev_mode.adjoint = 1; 
  ad_head(ad_outv,ad_inv,ad_vel,N);

#ifdef _CIVL
  for (i = 0; i < N; i++) global_ad_outv[i] = DERIV_val(ad_outv[i]); 
#endif

#ifdef DEBUG
  printf("Adjoint Output inv \n");
  double temp_adj;
  for (i = 0; i <N; i++) {  
    temp_adj = 0.0; 
    for (j = 0; j <ADIC_GRADVEC_LENGTH; j++) {
      temp_adj +=  DERIV_grad(ad_inv[i])[j];
    }
#ifdef _CIVL
    global_ad_inv_b[i] = temp_adj;
#endif
    printf("%d %lf \n", i, temp_adj); 
  } 
  printf("Adjoint Output vel \n");
  for (i = 0; i <N; i++) {  
    temp_adj = 0.0; 
    for (j = 0; j <ADIC_GRADVEC_LENGTH; j++) {
      temp_adj +=  DERIV_grad(ad_vel[i])[j];
    }
#ifdef _CIVL
    global_ad_vel_b[i] = temp_adj;
#endif
    printf("%d %lf \n", i, temp_adj); 
  }
#endif  
  ADIC_Finalize();
  free(ad_outv);
  free(ad_inv);
  free(ad_vel);
}

