#include <stdio.h>
#include "math.h"
#include "stdlib.h"


void an_head (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b, int n);
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

void head_b (double *outv, double *inv,double *vel, double *outv_b, double *inv_b,double *vel_b){
  int i;
  an_head(outv,inv, vel,outv_b,inv_b,vel_b, N);

  printf("Adjoint Output inv \n");
  double temp_adj;
  for (i = 0; i <N; i++) {  
    printf("%d %lf \n", i, inv_b[i]); 
  } 
  printf("Adjoint Output vel \n");
  for (i = 0; i <N; i++) {  
    printf("%d %lf \n", i, vel_b[i]); 
  }
}

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

