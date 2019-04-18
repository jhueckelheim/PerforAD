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
#endif
#ifdef ADIC_DENSE_REVERSE
#include "ad_types.hpp"
#endif
#ifdef ADIC_DENSE_SEED
#include "ad_types.h"
#endif
#ifdef ADIC_GRAD_LENGTH
#include "ad_grad_length_types.h"
#endif
#ifdef ADIC_SPARSE_NO_GRAD
#include "noderiv_sparslinc.h"
#endif
#ifdef ADIC_SPARSE
#include "sparslinc.h"
#endif
#define N 10
/*n)*/

void ad_head(DERIV_TYPE **outv,DERIV_TYPE **inv,DERIV_TYPE **vel,int n)
{
  int ad_Symbol_9;
  int ad_Symbol_8;
  int ad_Symbol_7;
  int ad_Symbol_6;
  int ad_Symbol_11;
  int ad_Symbol_10;
  int ad_Symbol_24;
  int ad_Symbol_23;
  int ad_Symbol_22;
  int ad_Symbol_21;
  DERIV_TYPE ad_prp_5;
  DERIV_TYPE ad_prp_4;
  DERIV_TYPE ad_prp_3;
  DERIV_TYPE ad_prp_2;
  DERIV_TYPE ad_prp_1;
  DERIV_TYPE ad_prp_0;
  DERIV_TYPE ad_prp_7;
  DERIV_TYPE ad_prp_6;
  double ad_Symbol_20;
  double ad_Symbol_19;
  double ad_Symbol_18;
  int ad_Symbol_17;
  int ad_Symbol_16;
  DERIV_TYPE ad_prp_8;
  int ad_Symbol_3;
  int ad_Symbol_2;
  int ad_Symbol_1;
  int ad_Symbol_0;
  int ad_Symbol_15;
  int ad_Symbol_14;
  int ad_Symbol_13;
  int ad_Symbol_12;
  double ad_lin_2;
  double ad_lin_1;
  double ad_lin_0;
  double ad_aux_1;
  double ad_aux_0;
  int ad_Symbol_5;
  int ad_Symbol_4;
  int ad_tyc_0;
  int ad_tyc_1;
  int ad_tyc_2;
  DERIV_TYPE ad_var_1;
  DERIV_TYPE *ad_var_0;
  int i;
  int j;
  ZeroDeriv(ad_prp_5);
  DERIV_val(ad_prp_5) = 0.00000;
  ZeroDeriv(ad_prp_4);
  DERIV_val(ad_prp_4) = 0.00000;
  ZeroDeriv(ad_prp_3);
  DERIV_val(ad_prp_3) = 0.00000;
  ZeroDeriv(ad_prp_2);
  DERIV_val(ad_prp_2) = 0.00000;
  ZeroDeriv(ad_prp_1);
  DERIV_val(ad_prp_1) = 0.00000;
  ZeroDeriv(ad_prp_0);
  DERIV_val(ad_prp_0) = 0.00000;
  ZeroDeriv(ad_prp_7);
  DERIV_val(ad_prp_7) = 0.00000;
  ZeroDeriv(ad_prp_6);
  DERIV_val(ad_prp_6) = 0.00000;
  ZeroDeriv(ad_prp_8);
  DERIV_val(ad_prp_8) = 0.00000;
  ZeroDeriv(ad_var_1);
  DERIV_val(ad_var_1) = 0.00000;
  if (our_rev_mode . plain == 1) {
    DERIV_TYPE ad_var_1;
    DERIV_TYPE *ad_var_0;
/* int i; int j; */
    int i;
    int j;
    for (j = 1; j <= n - 2; j = j + 1) {
      for (i = 1; i <= n - 2; i = i + 1) {
        ad_var_0 = &outv[i][j];
        DERIV_val(ad_var_1) = DERIV_val( *ad_var_0) + (- 4.00000 * DERIV_val(inv[i][j]) + DERIV_val(inv[i][j - 1]) + DERIV_val(inv[i][j + 1]) + DERIV_val(inv[i - 1][j]) + DERIV_val(inv[i + 1][j])) * DERIV_val(vel[i][j]);
        DERIV_val( *ad_var_0) = DERIV_val(ad_var_1);
      }
    }
  }
   else if (our_rev_mode . tape == 1) {
    DERIV_TYPE ad_var_1;
    DERIV_TYPE *ad_var_0;
/* int i; int j; */
    int i;
    int j;
    allocateBytesForSizeOf_SgTypeInt(&ad_tyc_2,1);
    FW_ADMM_Allocate(&j,ad_tyc_2,0);
    allocateBytesForSizeOf_SgTypeInt(&ad_tyc_1,1);
    FW_ADMM_Allocate(&i,ad_tyc_1,0);
    allocateBytesForSizeOf_SgTypeDouble(&ad_tyc_0,1);
    FW_ADMM_Allocate(&ad_var_1,ad_tyc_0,0);
    ad_Symbol_4 = 0;
    for (j = 1; j <= n - 2; j = j + 1) {
      ad_Symbol_5 = 0;
      for (i = 1; i <= n - 2; i = i + 1) {
        push_p_s0(ad_var_0);
        ad_var_0 = &outv[i][j];
        ad_aux_0 = - 4.00000;
        ad_aux_1 = ad_aux_0 * DERIV_val(inv[i][j]) + DERIV_val(inv[i][j - 1]) + DERIV_val(inv[i][j + 1]) + DERIV_val(inv[i - 1][j]) + DERIV_val(inv[i + 1][j]);
        ad_lin_0 = ad_aux_0;
        ad_lin_1 = DERIV_val(vel[i][j]);
        ad_lin_2 = ad_aux_1;
        DERIV_val(ad_var_1) = DERIV_val( *ad_var_0) + ad_aux_1 * DERIV_val(vel[i][j]);
        ad_Symbol_12 = j - 1;
        push_i_s0(ad_Symbol_12);
        ad_Symbol_13 = j + 1;
        push_i_s0(ad_Symbol_13);
        ad_Symbol_14 = i - 1;
        push_i_s0(ad_Symbol_14);
        ad_Symbol_15 = i + 1;
        push_i_s0(ad_Symbol_15);
        push_s0(ad_lin_0);
        push_s0(ad_lin_1);
        push_s0(ad_lin_2);
        push_i_s0(i);
        push_i_s0(j);
        DERIV_val( *ad_var_0) = DERIV_val(ad_var_1);
        ad_Symbol_5 = ad_Symbol_5 + 1;
      }
      push_i_s0(ad_Symbol_5);
      ad_Symbol_4 = ad_Symbol_4 + 1;
    }
    push_i_s0(ad_Symbol_4);
    push_p_s0(ad_var_0);
    FW_ADMM_Deallocate(&ad_var_1,0,0);
    FW_ADMM_Deallocate(&i,0,0);
    FW_ADMM_Deallocate(&j,0,0);
  }
   else if (our_rev_mode . adjoint == 1) {
    BW_ADMM_Deallocate(&j,0,0,0);
    BW_ADMM_Deallocate(&i,0,0,0);
    BW_ADMM_Deallocate(&ad_var_1,0,0,0);
    pop_p_s0(ad_var_0);
    ADMM_Rebase(&ad_var_0);
    pop_i_s0(ad_Symbol_0);
    for (ad_Symbol_1 = 1; ad_Symbol_1 <= ad_Symbol_0; ad_Symbol_1 = ad_Symbol_1 + 1) {
      pop_i_s0(ad_Symbol_2);
      for (ad_Symbol_3 = 1; ad_Symbol_3 <= ad_Symbol_2; ad_Symbol_3 = ad_Symbol_3 + 1) {
        IncDeriv(ad_prp_8, *ad_var_0);
        ZeroDeriv( *ad_var_0);
        IncDeriv(ad_var_1,ad_prp_8);
        ZeroDeriv(ad_prp_8);
        pop_i_s0(ad_Symbol_16);
        pop_i_s0(ad_Symbol_17);
        pop_s0(ad_Symbol_18);
        pop_s0(ad_Symbol_19);
        pop_s0(ad_Symbol_20);
        Saxpy(ad_Symbol_18,ad_var_1,ad_prp_6);
        Saxpy(ad_Symbol_19,ad_var_1,ad_prp_7);
        IncDeriv(ad_prp_0,ad_var_1);
        ZeroDeriv(ad_var_1);
        Saxpy(ad_Symbol_20,ad_prp_7,ad_prp_1);
        IncDeriv(ad_prp_2,ad_prp_7);
        IncDeriv(ad_prp_3,ad_prp_7);
        IncDeriv(ad_prp_4,ad_prp_7);
        IncDeriv(ad_prp_5,ad_prp_7);
        ZeroDeriv(ad_prp_7);
        IncDeriv(vel[ad_Symbol_17][ad_Symbol_16],ad_prp_6);
        ZeroDeriv(ad_prp_6);
        pop_i_s0(ad_Symbol_21);
        IncDeriv(inv[ad_Symbol_21][ad_Symbol_16],ad_prp_5);
        ZeroDeriv(ad_prp_5);
        pop_i_s0(ad_Symbol_22);
        IncDeriv(inv[ad_Symbol_22][ad_Symbol_16],ad_prp_4);
        ZeroDeriv(ad_prp_4);
        pop_i_s0(ad_Symbol_23);
        IncDeriv(inv[ad_Symbol_17][ad_Symbol_23],ad_prp_3);
        ZeroDeriv(ad_prp_3);
        pop_i_s0(ad_Symbol_24);
        IncDeriv(inv[ad_Symbol_17][ad_Symbol_24],ad_prp_2);
        ZeroDeriv(ad_prp_2);
        IncDeriv(inv[ad_Symbol_17][ad_Symbol_16],ad_prp_1);
        ZeroDeriv(ad_prp_1);
        IncDeriv( *ad_var_0,ad_prp_0);
        ZeroDeriv(ad_prp_0);
        pop_p_s0(ad_var_0);
        ADMM_Rebase(&ad_var_0);
      }
    }
    BW_ADMM_Allocate(&ad_var_1,0);
    BW_ADMM_Allocate(&i,0);
    BW_ADMM_Allocate(&j,0);
  }
}
