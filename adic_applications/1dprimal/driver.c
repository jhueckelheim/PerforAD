#include <stdio.h>

#include "ad_types.h"
#include "math.h"
#include "stdlib.h"

void head (double *outv, double *invec,double *vel, int n);
void ad_head(DERIV_TYPE *outv,DERIV_TYPE *invec,DERIV_TYPE *vel,int n);

#ifdef DEBUG
void an_head (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b,int n);
#endif

int main()
{
  double outv[N], inv[N], vel[N];
  DERIV_TYPE ad_outv[N], ad_inv[N], ad_vel[N];
  int i,j;
#ifdef DEBUG
  double invph[N], outvph[N];
  double h, dd_inv[N];
  double outv_b[N], inv_b[N], vel_b[N];

  for (i = 0; i < N; i++){
    outv[i] = 0.0;
    inv[i] = (i+1)*(i+1) * 0.5;
    invph[i] = (i+1)*(i+1) * 0.5;
    vel[i] = (i+1)*(i+1) * 0.5;
  }

  /* dd step size */
  h=.0001;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      outv[j] = 0.0;
      outvph[j] = 0.0;
    }
    invph[i] = inv[i] + h; 

    head(outv,inv, vel,N);
    head(outvph,invph, vel,N);
    dd_inv[i] = 0.0;
    for (j = 0; j < N; j++) {
      dd_inv[i] += (outvph[j] - outv[j])/h;
    }
    invph[i] = inv[i];
  } 
  for (i = 0; i < N; i++){
    outv[i] = 0.0;
    inv[i] = (i+1)*(i+1) * 0.5;
    vel[i] = (i+1)*(i+1) * 0.5;
  }
  head(outv,inv, vel,N);
  for (i = 0; i < N; i++){
    outv_b[i] = 1.0;
    inv_b[i] = 0.0;
    vel_b[i] = 0.0;
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
  ADIC_SetDepArray(ad_outv, N);
  ADIC_SetIndepArray(ad_inv, N);
  ADIC_SetIndepArray(ad_vel, N);
  ADIC_SetIndepDone();
  // Initialize the value of the independent variable ad_x
  for (i = 0; i < N; i++){
    inv[i] = (i+1)*(i+1) * 0.5;
    vel[i] = (i+1)*(i+1) * 0.5;
    DERIV_val(ad_outv[i]) =0.0;
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
    printf("%d\t%lf\t%lf \n", i, outv[i], DERIV_val(ad_outv[i]));
   } 
#endif
  our_rev_mode.tape = 0; 
  our_rev_mode.adjoint = 1; 
  ad_head(ad_outv,ad_inv,ad_vel,N);

#ifdef DEBUG
  printf("Adjoint Output inv \n");
  printf("i\tADIC\tPerforAD\tFD\n");
  double temp_adj;
  for (i = 0; i <N; i++) {  
    temp_adj = 0.0; 
    for (j = 0; j <ADIC_GRADVEC_LENGTH; j++) {
      temp_adj +=  DERIV_grad(ad_inv[i])[j];
    }
    printf("%d %lf\t%lf\t%lf \n", i,temp_adj, inv_b[i], dd_inv[i]); 
  } 
  printf("Adjoint Output vel \n");
  printf("i\tADIC\tPerforAD\n");
  for (i = 0; i <N; i++) {  
    temp_adj = 0.0; 
    for (j = 0; j <ADIC_GRADVEC_LENGTH; j++) {
      temp_adj +=  DERIV_grad(ad_vel[i])[j];
    }
    printf("%d %lf\t%lf \n", i, temp_adj, vel_b[i]); 
  }
#endif  
  ADIC_Finalize();

  return 0;
}
#ifdef DEBUG
void an_head (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b,int n){
int i=0;
i=0;
inv_b[i] += outv_b[i + 1]*vel[i + 1];
i=n - 2;
inv_b[i] += -2.0*outv_b[i]*vel[i];
vel_b[i]+= (-2.0*inv[i] + inv[i - 1] + inv[i + 1])*outv_b[i];
inv_b[i] += outv_b[i - 1]*vel[i - 1];
i=1;
inv_b[i] += outv_b[i + 1]*vel[i + 1];
inv_b[i] += -2.0*outv_b[i]*vel[i];
vel_b[i] += (-2.0*inv[i] + inv[i - 1] + inv[i + 1])*outv_b[i];
i=n - 1;
inv_b[i] += outv_b[i - 1]*vel[i - 1];
for ( i=2; i<=n - 3; i++ ) {
   inv_b[i] += outv_b[i + 1]*vel[i + 1];
   inv_b[i] += -2.0*outv_b[i]*vel[i];
   vel_b[i] += (-2.0*inv[i] + inv[i - 1] + inv[i + 1])*outv_b[i];
   inv_b[i] += outv_b[i - 1]*vel[i - 1];
}
}
#endif

