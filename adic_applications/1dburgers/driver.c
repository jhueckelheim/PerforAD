#include <stdio.h>

#include "ad_types.h"
#include "math.h"
#include "stdlib.h"


void head (double *u, double *u_1, double C, double D,int n);
void ad_head(DERIV_TYPE *u,DERIV_TYPE *u_1,DERIV_TYPE *C,DERIV_TYPE *D,int n);

#ifdef DEBUG
void an_head (double *u, double *u_1, double C, double D, double *u_b, double *u_1_b, int n);
#endif

int main()
{
  double u[N], u_1[N];
  DERIV_TYPE ad_u[N], ad_u_1[N], ad_C, ad_D;
  int i,j;
#ifdef DEBUG
  double u_1ph[N], uph[N];
  double h, dd_u_1[N];
  double u_b[N], u_1_b[N];
  double C, D;

  C=0.5;
  D=0.5;
  for (i = 0; i < N; i++){
    u[i] = 0.0;
    u_1[i] = (i+1)*(i+1) * 0.5;
    u_1ph[i] = (i+1)*(i+1) * 0.5;
  }

  /* dd step size */
  h=.0001;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      u[j] = 0.0;
      uph[j] = 0.0;
    }
    u_1ph[i] = u_1[i] + h; 

    head(u,u_1,C,D,N);
    head(uph,u_1ph,C,D,N);
    dd_u_1[i] = 0.0;
    for (j = 0; j < N; j++) {
      dd_u_1[i] += (uph[j] - u[j])/h;
    }
    u_1ph[i] = u_1[i];
  } 
  for (i = 0; i < N; i++){
    u[i] = 0.0;
    u_1[i] = (i+1)*(i+1) * 0.5;
  }
  head(u,u_1,C,D,N);
  for (i = 0; i < N; i++){
    u_b[i] = 1.0;
    u_1_b[i] = 0.0;
  }
  an_head(u,u_1,C,D,u_b,u_1_b, N);
#endif
  /* ----------------------------------------------------------
   * Compute dy/dx using AD 
   */
  ADIC_Init();
  __ADIC_TapeInit();
  
  
  // Set indpendent variables 
  ADIC_SetReverseMode();
  ADIC_SetDepArray(ad_u, N);
  ADIC_SetIndepArray(ad_u_1, N);
  ADIC_SetIndepDone();
  // Initialize the value of the independent variable ad_x
  for (i = 0; i < N; i++){
    u_1[i] = (i+1)*(i+1) * 0.5;
    DERIV_val(ad_u[i]) =0.0;
    DERIV_val(ad_u_1[i]) =u_1[i];
  }
  ZeroDeriv(ad_C);
  ZeroDeriv(ad_D);
  DERIV_val(ad_C) = 0.5;
  DERIV_val(ad_D) = 0.5;
   
  
  // Invoke AD function 
  our_rev_mode.tape = 1; 
  our_rev_mode.adjoint = 0; 
  ad_head(ad_u,ad_u_1,&ad_C,&ad_D,N);
#ifdef DEBUG
  printf("Primal Output\n");
  printf("i\tOrig\tADIC\n");
  for (i = 0; i < N; i++){
    printf("%d\t%lf\t%lf \n", i, u[i], DERIV_val(ad_u[i]));
   } 
#endif
  our_rev_mode.tape = 0; 
  our_rev_mode.adjoint = 1; 
  ad_head(ad_u,ad_u_1,&ad_C,&ad_D,N);

#ifdef DEBUG
  printf("Adjoint Output u_1 \n");
  printf("i\tADIC\tPerforAD\tFD\n");
  double temp_adj;
  for (i = 0; i <N; i++) {  
    temp_adj = 0.0; 
    for (j = 0; j <ADIC_GRADVEC_LENGTH; j++) {
      temp_adj +=  DERIV_grad(ad_u_1[i])[j];
    }
    printf("%d %lf\t%lf\t%lf \n", i,temp_adj, u_1_b[i], dd_u_1[i]); 
  } 
#endif  
  ADIC_Finalize();

  return 0;
}
#ifdef DEBUG
void an_head (double *u, double *u_1, double C, double D, double *u_b, double *u_1_b, int n){
int i=0;
i=0;
u_1_b[i] += (C*((0> u_1[i + 1])?0: u_1[i + 1]) + D)*u_b[i + 1];
i=n - 2;
u_1_b[i] += (-C*((-u_1[i] + u_1[i + 1])*((-u_1[i]>=0)?1.0:0.0) + (u_1[i] - u_1[i - 1])*((u_1[i]>=0)?1.0:0.0) + ((0> u_1[i])?0: u_1[i]) - ((0< u_1[i])?0: u_1[i])) - 2.0*D + 1)*u_b[i];
u_1_b[i] += (-C*((0< u_1[i - 1])?0: u_1[i - 1]) + D)*u_b[i - 1];
i=1;
u_1_b[i] += (C*((0> u_1[i + 1])?0: u_1[i + 1]) + D)*u_b[i + 1];
u_1_b[i] += (-C*((-u_1[i] + u_1[i + 1])*((-u_1[i]>=0)?1.0:0.0) + (u_1[i] - u_1[i - 1])*((u_1[i]>=0)?1.0:0.0) + ((0> u_1[i])?0: u_1[i]) - ((0< u_1[i])?0: u_1[i])) - 2.0*D + 1)*u_b[i];
i=n - 1;
u_1_b[i] += (-C*((0< u_1[i - 1])?0: u_1[i - 1]) + D)*u_b[i - 1];
for ( i=2; i<=n - 3; i++ ) {
    u_1_b[i] += (C*((0> u_1[i + 1])?0: u_1[i + 1]) + D)*u_b[i + 1];
    u_1_b[i] += (-C*((-u_1[i] + u_1[i + 1])*((-u_1[i]>=0)?1.0:0.0) + (u_1[i] - u_1[i - 1])*((u_1[i]>=0)?1.0:0.0) + ((0> u_1[i])?0: u_1[i]) - ((0< u_1[i])?0: u_1[i])) - 2.0*D + 1)*u_b[i];
    u_1_b[i] += (-C*((0< u_1[i - 1])?0: u_1[i - 1]) + D)*u_b[i - 1];
}

}
#endif

