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
/*n)*/

void ad_head(DERIV_TYPE ***u,DERIV_TYPE ***u_1,DERIV_TYPE ***u_2,DERIV_TYPE ***c,DERIV_TYPE *D,int n)
{
  int ad_Symbol_14;
  int ad_Symbol_13;
  int ad_Symbol_12;
  int ad_Symbol_11;
  int ad_Symbol_10;
  int ad_Symbol_9;
  int ad_Symbol_17;
  int ad_Symbol_16;
  int ad_Symbol_15;
  int ad_Symbol_36;
  int ad_Symbol_35;
  int ad_Symbol_34;
  int ad_Symbol_33;
  int ad_Symbol_32;
  int ad_Symbol_31;
  DERIV_TYPE ad_prp_8;
  DERIV_TYPE ad_prp_7;
  DERIV_TYPE ad_prp_6;
  DERIV_TYPE ad_prp_5;
  DERIV_TYPE ad_prp_4;
  DERIV_TYPE ad_prp_3;
  DERIV_TYPE ad_prp_2;
  DERIV_TYPE ad_prp_0;
  DERIV_TYPE ad_prp_11;
  DERIV_TYPE ad_prp_10;
  DERIV_TYPE ad_prp_9;
  DERIV_TYPE ad_prp_1;
  DERIV_TYPE ad_prp_12;
  int ad_Symbol_30;
  double ad_Symbol_29;
  double ad_Symbol_28;
  double ad_Symbol_27;
  int ad_Symbol_26;
  int ad_Symbol_25;
  int ad_Symbol_24;
  DERIV_TYPE ad_prp_13;
  int ad_Symbol_5;
  int ad_Symbol_4;
  int ad_Symbol_3;
  int ad_Symbol_2;
  int ad_Symbol_1;
  int ad_Symbol_0;
  int ad_Symbol_23;
  int ad_Symbol_22;
  int ad_Symbol_21;
  int ad_Symbol_20;
  int ad_Symbol_19;
  int ad_Symbol_18;
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
  int ad_Symbol_8;
  int ad_Symbol_7;
  int ad_Symbol_6;
  int ad_tyc_0;
  int ad_tyc_1;
  int ad_tyc_2;
  int ad_tyc_3;
  DERIV_TYPE ad_var_1;
  DERIV_TYPE *ad_var_0;
  int i;
  int j;
  int k;
  ZeroDeriv(ad_prp_8);
  DERIV_val(ad_prp_8) = 0.00000;
  ZeroDeriv(ad_prp_7);
  DERIV_val(ad_prp_7) = 0.00000;
  ZeroDeriv(ad_prp_6);
  DERIV_val(ad_prp_6) = 0.00000;
  ZeroDeriv(ad_prp_5);
  DERIV_val(ad_prp_5) = 0.00000;
  ZeroDeriv(ad_prp_4);
  DERIV_val(ad_prp_4) = 0.00000;
  ZeroDeriv(ad_prp_3);
  DERIV_val(ad_prp_3) = 0.00000;
  ZeroDeriv(ad_prp_2);
  DERIV_val(ad_prp_2) = 0.00000;
  ZeroDeriv(ad_prp_0);
  DERIV_val(ad_prp_0) = 0.00000;
  ZeroDeriv(ad_prp_11);
  DERIV_val(ad_prp_11) = 0.00000;
  ZeroDeriv(ad_prp_10);
  DERIV_val(ad_prp_10) = 0.00000;
  ZeroDeriv(ad_prp_9);
  DERIV_val(ad_prp_9) = 0.00000;
  ZeroDeriv(ad_prp_1);
  DERIV_val(ad_prp_1) = 0.00000;
  ZeroDeriv(ad_prp_12);
  DERIV_val(ad_prp_12) = 0.00000;
  ZeroDeriv(ad_prp_13);
  DERIV_val(ad_prp_13) = 0.00000;
  ZeroDeriv(ad_var_1);
  DERIV_val(ad_var_1) = 0.00000;
  if (our_rev_mode . plain == 1) {
    DERIV_TYPE ad_var_1;
    DERIV_TYPE *ad_var_0;
/* int i; int j; int k; */
    int i;
    int j;
    int k;
    for (k = 1; k <= n - 2; k = k + 1) {
      for (j = 1; j <= n - 2; j = j + 1) {
        for (i = 1; i <= n - 2; i = i + 1) {
          ad_var_0 = &u[i][j][k];
          DERIV_val(ad_var_1) = DERIV_val( *ad_var_0) + (DERIV_val( *D) * (- 6 * DERIV_val(u_1[i][j][k]) + DERIV_val(u_1[i][j][k - 1]) + DERIV_val(u_1[i][j][k + 1]) + DERIV_val(u_1[i][j - 1][k]) + DERIV_val(u_1[i][j + 1][k]) + DERIV_val(u_1[i - 1][j][k]) + DERIV_val(u_1[i + 1][j][k])) * DERIV_val(c[i][j][k]) + 2.00000 * DERIV_val(u_1[i][j][k]) - DERIV_val(u_2[i][j][k]));
          DERIV_val( *ad_var_0) = DERIV_val(ad_var_1);
        }
      }
    }
  }
   else if (our_rev_mode . tape == 1) {
    DERIV_TYPE ad_var_1;
    DERIV_TYPE *ad_var_0;
/* int i; int j; int k; */
    int i;
    int j;
    int k;
    allocateBytesForSizeOf_SgTypeInt(&ad_tyc_3,1);
    FW_ADMM_Allocate(&k,ad_tyc_3,0);
    allocateBytesForSizeOf_SgTypeInt(&ad_tyc_2,1);
    FW_ADMM_Allocate(&j,ad_tyc_2,0);
    allocateBytesForSizeOf_SgTypeInt(&ad_tyc_1,1);
    FW_ADMM_Allocate(&i,ad_tyc_1,0);
    allocateBytesForSizeOf_SgTypeDouble(&ad_tyc_0,1);
    FW_ADMM_Allocate(&ad_var_1,ad_tyc_0,0);
    ad_Symbol_6 = 0;
    for (k = 1; k <= n - 2; k = k + 1) {
      ad_Symbol_7 = 0;
      for (j = 1; j <= n - 2; j = j + 1) {
        ad_Symbol_8 = 0;
        for (i = 1; i <= n - 2; i = i + 1) {
          push_p_s0(ad_var_0);
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
          ad_Symbol_18 = k - 1;
          push_i_s0(ad_Symbol_18);
          ad_Symbol_19 = k + 1;
          push_i_s0(ad_Symbol_19);
          ad_Symbol_20 = j - 1;
          push_i_s0(ad_Symbol_20);
          ad_Symbol_21 = j + 1;
          push_i_s0(ad_Symbol_21);
          ad_Symbol_22 = i - 1;
          push_i_s0(ad_Symbol_22);
          ad_Symbol_23 = i + 1;
          push_i_s0(ad_Symbol_23);
          push_s0(ad_lin_0);
          push_s0(ad_lin_4);
          push_s0(ad_acc_0);
          push_s0(ad_acc_1);
          push_i_s0(i);
          push_i_s0(j);
          push_i_s0(k);
          DERIV_val( *ad_var_0) = DERIV_val(ad_var_1);
          ad_Symbol_8 = ad_Symbol_8 + 1;
        }
        push_i_s0(ad_Symbol_8);
        ad_Symbol_7 = ad_Symbol_7 + 1;
      }
      push_i_s0(ad_Symbol_7);
      ad_Symbol_6 = ad_Symbol_6 + 1;
    }
    push_i_s0(ad_Symbol_6);
    push_p_s0(ad_var_0);
    FW_ADMM_Deallocate(&ad_var_1,0,0);
    FW_ADMM_Deallocate(&i,0,0);
    FW_ADMM_Deallocate(&j,0,0);
    FW_ADMM_Deallocate(&k,0,0);
  }
   else if (our_rev_mode . adjoint == 1) {
    BW_ADMM_Deallocate(&k,0,0,0);
    BW_ADMM_Deallocate(&j,0,0,0);
    BW_ADMM_Deallocate(&i,0,0,0);
    BW_ADMM_Deallocate(&ad_var_1,0,0,0);
    pop_p_s0(ad_var_0);
    ADMM_Rebase(&ad_var_0);
    pop_i_s0(ad_Symbol_0);
    for (ad_Symbol_1 = 1; ad_Symbol_1 <= ad_Symbol_0; ad_Symbol_1 = ad_Symbol_1 + 1) {
      pop_i_s0(ad_Symbol_2);
      for (ad_Symbol_3 = 1; ad_Symbol_3 <= ad_Symbol_2; ad_Symbol_3 = ad_Symbol_3 + 1) {
        pop_i_s0(ad_Symbol_4);
        for (ad_Symbol_5 = 1; ad_Symbol_5 <= ad_Symbol_4; ad_Symbol_5 = ad_Symbol_5 + 1) {
          IncDeriv(ad_prp_13, *ad_var_0);
          ZeroDeriv( *ad_var_0);
          IncDeriv(ad_var_1,ad_prp_13);
          ZeroDeriv(ad_prp_13);
          pop_i_s0(ad_Symbol_24);
          pop_i_s0(ad_Symbol_25);
          pop_i_s0(ad_Symbol_26);
          pop_s0(ad_Symbol_27);
          pop_s0(ad_Symbol_28);
          pop_s0(ad_Symbol_29);
          pop_s0(ad_Symbol_30);
          Saxpy(ad_Symbol_27,ad_var_1,ad_prp_12);
          Saxpy(ad_Symbol_28,ad_var_1,ad_prp_1);
          Saxpy(ad_Symbol_29,ad_var_1,ad_prp_9);
          Saxpy(2.00000,ad_var_1,ad_prp_10);
          DecDeriv(ad_prp_11,ad_var_1);
          IncDeriv(ad_prp_0,ad_var_1);
          ZeroDeriv(ad_var_1);
          Saxpy(ad_Symbol_30,ad_prp_12,ad_prp_2);
          IncDeriv(ad_prp_3,ad_prp_12);
          IncDeriv(ad_prp_4,ad_prp_12);
          IncDeriv(ad_prp_5,ad_prp_12);
          IncDeriv(ad_prp_6,ad_prp_12);
          IncDeriv(ad_prp_7,ad_prp_12);
          IncDeriv(ad_prp_8,ad_prp_12);
          ZeroDeriv(ad_prp_12);
          IncDeriv(u_2[ad_Symbol_26][ad_Symbol_25][ad_Symbol_24],ad_prp_11);
          ZeroDeriv(ad_prp_11);
          IncDeriv(u_1[ad_Symbol_26][ad_Symbol_25][ad_Symbol_24],ad_prp_10);
          ZeroDeriv(ad_prp_10);
          IncDeriv(c[ad_Symbol_26][ad_Symbol_25][ad_Symbol_24],ad_prp_9);
          ZeroDeriv(ad_prp_9);
          pop_i_s0(ad_Symbol_31);
          IncDeriv(u_1[ad_Symbol_31][ad_Symbol_25][ad_Symbol_24],ad_prp_8);
          ZeroDeriv(ad_prp_8);
          pop_i_s0(ad_Symbol_32);
          IncDeriv(u_1[ad_Symbol_32][ad_Symbol_25][ad_Symbol_24],ad_prp_7);
          ZeroDeriv(ad_prp_7);
          pop_i_s0(ad_Symbol_33);
          IncDeriv(u_1[ad_Symbol_26][ad_Symbol_33][ad_Symbol_24],ad_prp_6);
          ZeroDeriv(ad_prp_6);
          pop_i_s0(ad_Symbol_34);
          IncDeriv(u_1[ad_Symbol_26][ad_Symbol_34][ad_Symbol_24],ad_prp_5);
          ZeroDeriv(ad_prp_5);
          pop_i_s0(ad_Symbol_35);
          IncDeriv(u_1[ad_Symbol_26][ad_Symbol_25][ad_Symbol_35],ad_prp_4);
          ZeroDeriv(ad_prp_4);
          pop_i_s0(ad_Symbol_36);
          IncDeriv(u_1[ad_Symbol_26][ad_Symbol_25][ad_Symbol_36],ad_prp_3);
          ZeroDeriv(ad_prp_3);
          IncDeriv(u_1[ad_Symbol_26][ad_Symbol_25][ad_Symbol_24],ad_prp_2);
          ZeroDeriv(ad_prp_2);
          IncDeriv( *D,ad_prp_1);
          ZeroDeriv(ad_prp_1);
          IncDeriv( *ad_var_0,ad_prp_0);
          ZeroDeriv(ad_prp_0);
          pop_p_s0(ad_var_0);
          ADMM_Rebase(&ad_var_0);
        }
      }
    }
    BW_ADMM_Allocate(&ad_var_1,0);
    BW_ADMM_Allocate(&i,0);
    BW_ADMM_Allocate(&j,0);
    BW_ADMM_Allocate(&k,0);
  }
}
