/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 04/17/19 21:31:40 by                  */
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
#define N 10
/*n)*/

void ad_head(DERIV_TYPE **outv,DERIV_TYPE **inv,DERIV_TYPE **vel,int n)
{
  DERIV_TYPE ad_var_1;
  DERIV_TYPE *ad_var_0;
/* int i; int j; */
  int i;
  int j;
/*for(j = 1;j <= n - 2;j++)*/
  for (j = 1, j = 1; j <= n - 2; j = j + 1) {
/*for(i = 1;i <= n - 2;i++)*/
    for (i = 1, i = 1; i <= n - 2; i = i + 1) {
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
      double ad_lin_2;
      double ad_lin_1;
      double ad_lin_0;
      double ad_aux_1;
      double ad_aux_0;
      ad_var_0 = &outv[i][j];
      ad_aux_0 = - 4.00000;
      ad_aux_1 = ad_aux_0 * DERIV_val(inv[i][j]) + DERIV_val(inv[i][j - 1]) + DERIV_val(inv[i][j + 1]) + DERIV_val(inv[i - 1][j]) + DERIV_val(inv[i + 1][j]);
      ad_lin_0 = ad_aux_0;
      ad_lin_1 = DERIV_val(vel[i][j]);
      ad_lin_2 = ad_aux_1;
      DERIV_val(ad_var_1) = DERIV_val( *ad_var_0) + ad_aux_1 * DERIV_val(vel[i][j]);
      ADIC_SetDeriv( *ad_var_0,ad_prp_0);
      ADIC_SetDeriv(inv[i][j],ad_prp_1);
      ADIC_SetDeriv(inv[i][j - 1],ad_prp_2);
      ADIC_SetDeriv(inv[i][j + 1],ad_prp_3);
      ADIC_SetDeriv(inv[i - 1][j],ad_prp_4);
      ADIC_SetDeriv(inv[i + 1][j],ad_prp_5);
      ADIC_SetDeriv(vel[i][j],ad_prp_6);
      ADIC_Sax_5(1,ad_prp_5,1,ad_prp_4,1,ad_prp_3,1,ad_prp_2,ad_lin_0,ad_prp_1,ad_prp_7);
      ADIC_Sax_3(1,ad_prp_0,ad_lin_1,ad_prp_7,ad_lin_2,ad_prp_6,ad_var_1);
      DERIV_val( *ad_var_0) = DERIV_val(ad_var_1);
      ADIC_SetDeriv(ad_var_1,ad_prp_8);
      ADIC_Sax_1(1,ad_prp_8, *ad_var_0);
    }
{
    }
  }
{
  }
  return ;
}
