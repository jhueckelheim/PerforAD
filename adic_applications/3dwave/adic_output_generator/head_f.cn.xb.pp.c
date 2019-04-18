/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 04/18/19 00:03:13 by                  */
/*   ADIC version /*                                                                    */
/*   ADIC was prepared as an account of work sponsored by an          */
/*   agency of the United States Government and the University of     */
/*   Chicago.  NEITHER THE AUTHOR(S), THE UNITED STATES GOVERNMENT    */
/*   NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, INCLUDING */
/*   ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS  */
/*   OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR */
/*   THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR  */
/*   PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE */
/*   PRIVATELY OWNED RIGHTS.                                          */
/*                                                                    */
/**********************************************************************/
#ifdef ADIC_DENSE
#include "ad_types.h"
#include "ad_grad_f_saxpy-n_dense.h"
#endif
#ifdef ADIC_DENSE_REVERSE
#include "ad_types.hpp"
#include "ad_grad_f_saxpy-n_dense.h"
#endif
#ifdef ADIC_DENSE_SEED
#include "ad_types.h"
#include "ad_grad_f_saxpy-n_dense.h"
#endif
#ifdef ADIC_GRAD_LENGTH
#include "ad_grad_length_types.h"
#include "ad_grad_f_saxpy-n_dense.h"
#endif
#ifdef ADIC_SPARSE_NO_GRAD
#include "noderiv_sparslinc.h"
#include "ad_grad_saxpy-n_sparse.h"
#endif
#ifdef ADIC_SPARSE
#include "sparslinc.h"
#include "ad_grad_saxpy-n_sparse.h"
#endif
/*n)*/

void ad_head(DERIV_TYPE ***u,DERIV_TYPE ***u_1,DERIV_TYPE ***u_2,DERIV_TYPE ***c,DERIV_TYPE *D,int n)
{
  DERIV_TYPE ad_var_1;
  DERIV_TYPE *ad_var_0;
/* int i; int j; int k; */
  int i;
  int j;
  int k;
/*for(k = 1;k <= n - 2;k++)*/
  for (k = 1, k = 1; k <= n - 2; k = k + 1) {
/*for(j = 1;j <= n - 2;j++)*/
    for (j = 1, j = 1; j <= n - 2; j = j + 1) {
/*for(i = 1;i <= n - 2;i++)*/
      for (i = 1, i = 1; i <= n - 2; i = i + 1) {
        DERIV_TYPE ad_prp_13;
        ADIC_Initialize(&ad_prp_13);
        DERIV_TYPE ad_prp_12;
        ADIC_Initialize(&ad_prp_12);
        DERIV_TYPE ad_prp_11;
        ADIC_Initialize(&ad_prp_11);
        DERIV_TYPE ad_prp_10;
        ADIC_Initialize(&ad_prp_10);
        DERIV_TYPE ad_prp_9;
        ADIC_Initialize(&ad_prp_9);
        DERIV_TYPE ad_prp_8;
        ADIC_Initialize(&ad_prp_8);
        DERIV_TYPE ad_prp_7;
        ADIC_Initialize(&ad_prp_7);
        DERIV_TYPE ad_prp_6;
        ADIC_Initialize(&ad_prp_6);
        DERIV_TYPE ad_prp_5;
        ADIC_Initialize(&ad_prp_5);
        DERIV_TYPE ad_prp_4;
        ADIC_Initialize(&ad_prp_4);
        DERIV_TYPE ad_prp_3;
        ADIC_Initialize(&ad_prp_3);
        DERIV_TYPE ad_prp_2;
        ADIC_Initialize(&ad_prp_2);
        DERIV_TYPE ad_prp_1;
        ADIC_Initialize(&ad_prp_1);
        DERIV_TYPE ad_prp_0;
        ADIC_Initialize(&ad_prp_0);
        double ad_acc_1;
        double ad_acc_0;
        double ad_lin_4;
        double ad_lin_3;
        double ad_lin_2;
        int ad_lin_0;
        double ad_lin_1;
        double ad_aux_2;
        double ad_aux_1;
        int ad_aux_0;
        ad_var_0 = &u[i][j][k];
        ad_aux_0 = - 6;
        ad_aux_1 = ad_aux_0 * DERIV_val(u_1[i][j][k]) + DERIV_val(u_1[i][j][k - 1]) + DERIV_val(u_1[i][j][k + 1]) + DERIV_val(u_1[i][j - 1][k]) + DERIV_val(u_1[i][j + 1][k]) + DERIV_val(u_1[i - 1][j][k]) + DERIV_val(u_1[i + 1][j][k]);
        ad_aux_2 = DERIV_val( *D) * ad_aux_1;
        ad_lin_1 = ad_aux_1;
        ad_lin_0 = ad_aux_0;
        ad_lin_2 = DERIV_val( *D);
        ad_lin_3 = DERIV_val(c[i][j][k]);
        ad_lin_4 = ad_aux_2;
        DERIV_val(ad_var_1) = DERIV_val( *ad_var_0) + (ad_aux_2 * DERIV_val(c[i][j][k]) + 2.00000 * DERIV_val(u_1[i][j][k]) - DERIV_val(u_2[i][j][k]));
        ad_acc_0 = ad_lin_1 * ad_lin_3;
        ad_acc_1 = ad_lin_2 * ad_lin_3;
        ADIC_SetDeriv( *ad_var_0,ad_prp_0);
        ADIC_SetDeriv( *D,ad_prp_1);
        ADIC_SetDeriv(u_1[i][j][k],ad_prp_2);
        ADIC_SetDeriv(u_1[i][j][k - 1],ad_prp_3);
        ADIC_SetDeriv(u_1[i][j][k + 1],ad_prp_4);
        ADIC_SetDeriv(u_1[i][j - 1][k],ad_prp_5);
        ADIC_SetDeriv(u_1[i][j + 1][k],ad_prp_6);
        ADIC_SetDeriv(u_1[i - 1][j][k],ad_prp_7);
        ADIC_SetDeriv(u_1[i + 1][j][k],ad_prp_8);
        ADIC_SetDeriv(c[i][j][k],ad_prp_9);
        ADIC_SetDeriv(u_1[i][j][k],ad_prp_10);
        ADIC_SetDeriv(u_2[i][j][k],ad_prp_11);
        ADIC_Sax_7(1,ad_prp_8,1,ad_prp_7,1,ad_prp_6,1,ad_prp_5,1,ad_prp_4,1,ad_prp_3,ad_lin_0,ad_prp_2,ad_prp_12);
        ADIC_Sax_6(1,ad_prp_0,-1,ad_prp_11,2.00000,ad_prp_10,ad_lin_4,ad_prp_9,ad_acc_0,ad_prp_1,ad_acc_1,ad_prp_12,ad_var_1);
        DERIV_val( *ad_var_0) = DERIV_val(ad_var_1);
        ADIC_SetDeriv(ad_var_1,ad_prp_13);
        ADIC_Sax_1(1,ad_prp_13, *ad_var_0);
      }
{
      }
    }
{
    }
  }
{
  }
  return ;
}
