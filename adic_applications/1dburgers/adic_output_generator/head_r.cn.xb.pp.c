/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 04/17/19 23:26:05 by                  */
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

void ad_head(DERIV_TYPE *u,DERIV_TYPE *u_1,DERIV_TYPE *C,DERIV_TYPE *D,int n)
{
  int ad_Symbol_12;
  int ad_Symbol_11;
  int ad_Symbol_10;
  int ad_Symbol_9;
  int ad_Symbol_17;
  int ad_Symbol_16;
  int ad_Symbol_15;
  int ad_Symbol_14;
  int ad_Symbol_13;
  DERIV_TYPE ad_prp_19;
  int ad_Symbol_37;
  int ad_Symbol_3;
  DERIV_TYPE ad_prp_18;
  int ad_Symbol_36;
  int ad_Symbol_2;
  int ad_Symbol_35;
  int ad_Symbol_34;
  int ad_Symbol_33;
  int ad_Symbol_32;
  DERIV_TYPE ad_prp_3;
  DERIV_TYPE ad_prp_2;
  DERIV_TYPE ad_prp_5;
  DERIV_TYPE ad_prp_6;
  DERIV_TYPE ad_prp_13;
  DERIV_TYPE ad_prp_4;
  DERIV_TYPE ad_prp_14;
  DERIV_TYPE ad_prp_7;
  DERIV_TYPE ad_prp_11;
  DERIV_TYPE ad_prp_10;
  DERIV_TYPE ad_prp_9;
  DERIV_TYPE ad_prp_0;
  DERIV_TYPE ad_prp_12;
  DERIV_TYPE ad_prp_15;
  DERIV_TYPE ad_prp_1;
  DERIV_TYPE ad_prp_8;
  DERIV_TYPE ad_prp_16;
  double ad_Symbol_31;
  double ad_Symbol_30;
  double ad_Symbol_29;
  double ad_Symbol_28;
  double ad_Symbol_27;
  double ad_Symbol_26;
  double ad_Symbol_25;
  double ad_Symbol_24;
  double ad_Symbol_23;
  int ad_Symbol_22;
  DERIV_TYPE ad_prp_17;
  int ad_Symbol_1;
  int ad_Symbol_0;
  int ad_Symbol_21;
  int ad_Symbol_20;
  int ad_Symbol_19;
  int ad_Symbol_18;
  double ad_acc_0;
  double ad_lin_8;
  double ad_lin_6;
  double ad_lin_7;
  double ad_lin_5;
  double ad_lin_3;
  double ad_lin_2;
  double ad_lin_1;
  double ad_lin_0;
  double ad_lin_4;
  double ad_aux_5;
  double ad_aux_4;
  double ad_aux_3;
  double ad_aux_1;
  double ad_aux_0;
  double ad_aux_2;
  int ad_Symbol_8;
  int ad_Symbol_7;
  int ad_Symbol_6;
  int ad_Symbol_5;
  int ad_Symbol_4;
  int ad_tyc_0;
  int ad_tyc_1;
  int ad_tyc_2;
  int ad_tyc_3;
  DERIV_TYPE ad_var_3;
  DERIV_TYPE *ad_var_2;
  DERIV_TYPE ad_var_1;
  DERIV_TYPE ad_var_0;
  int i;
  ZeroDeriv(ad_prp_19);
  DERIV_val(ad_prp_19) = 0.00000;
  ZeroDeriv(ad_prp_18);
  DERIV_val(ad_prp_18) = 0.00000;
  ZeroDeriv(ad_prp_3);
  DERIV_val(ad_prp_3) = 0.00000;
  ZeroDeriv(ad_prp_2);
  DERIV_val(ad_prp_2) = 0.00000;
  ZeroDeriv(ad_prp_5);
  DERIV_val(ad_prp_5) = 0.00000;
  ZeroDeriv(ad_prp_6);
  DERIV_val(ad_prp_6) = 0.00000;
  ZeroDeriv(ad_prp_13);
  DERIV_val(ad_prp_13) = 0.00000;
  ZeroDeriv(ad_prp_4);
  DERIV_val(ad_prp_4) = 0.00000;
  ZeroDeriv(ad_prp_14);
  DERIV_val(ad_prp_14) = 0.00000;
  ZeroDeriv(ad_prp_7);
  DERIV_val(ad_prp_7) = 0.00000;
  ZeroDeriv(ad_prp_11);
  DERIV_val(ad_prp_11) = 0.00000;
  ZeroDeriv(ad_prp_10);
  DERIV_val(ad_prp_10) = 0.00000;
  ZeroDeriv(ad_prp_9);
  DERIV_val(ad_prp_9) = 0.00000;
  ZeroDeriv(ad_prp_0);
  DERIV_val(ad_prp_0) = 0.00000;
  ZeroDeriv(ad_prp_12);
  DERIV_val(ad_prp_12) = 0.00000;
  ZeroDeriv(ad_prp_15);
  DERIV_val(ad_prp_15) = 0.00000;
  ZeroDeriv(ad_prp_1);
  DERIV_val(ad_prp_1) = 0.00000;
  ZeroDeriv(ad_prp_8);
  DERIV_val(ad_prp_8) = 0.00000;
  ZeroDeriv(ad_prp_16);
  DERIV_val(ad_prp_16) = 0.00000;
  ZeroDeriv(ad_prp_17);
  DERIV_val(ad_prp_17) = 0.00000;
  ZeroDeriv(ad_var_3);
  DERIV_val(ad_var_3) = 0.00000;
  ZeroDeriv(ad_var_1);
  DERIV_val(ad_var_1) = 0.00000;
  ZeroDeriv(ad_var_0);
  DERIV_val(ad_var_0) = 0.00000;
  if (our_rev_mode . plain == 1) {
    DERIV_TYPE ad_var_3;
    DERIV_TYPE *ad_var_2;
    DERIV_TYPE ad_var_1;
    DERIV_TYPE ad_var_0;
/* int i; */
    int i;
    for (i = 1; i <= n - 2; i = i + 1) {
      if (0 < DERIV_val(u_1[i])) {
        DERIV_val(ad_var_0) = 0;
      }
       else {
        DERIV_val(ad_var_0) = DERIV_val(u_1[i]);
      }
      if (0 > DERIV_val(u_1[i])) {
        DERIV_val(ad_var_1) = 0;
      }
       else {
        DERIV_val(ad_var_1) = DERIV_val(u_1[i]);
      }
      ad_var_2 = &u[i];
      DERIV_val(ad_var_3) = DERIV_val( *ad_var_2) + (-DERIV_val( *C) * ((-DERIV_val(u_1[i]) + DERIV_val(u_1[i + 1])) * DERIV_val(ad_var_0) + (DERIV_val(u_1[i]) - DERIV_val(u_1[i - 1])) * DERIV_val(ad_var_1)) + DERIV_val( *D) * (- 2.00000 * DERIV_val(u_1[i]) + DERIV_val(u_1[i - 1]) + DERIV_val(u_1[i + 1])) + DERIV_val(u_1[i]));
      DERIV_val( *ad_var_2) = DERIV_val(ad_var_3);
    }
  }
   else if (our_rev_mode . tape == 1) {
    DERIV_TYPE ad_var_3;
    DERIV_TYPE *ad_var_2;
    DERIV_TYPE ad_var_1;
    DERIV_TYPE ad_var_0;
/* int i; */
    int i;
    allocateBytesForSizeOf_SgTypeInt(&ad_tyc_3,1);
    FW_ADMM_Allocate(&i,ad_tyc_3,0);
    allocateBytesForSizeOf_SgTypeDouble(&ad_tyc_2,1);
    FW_ADMM_Allocate(&ad_var_3,ad_tyc_2,0);
    allocateBytesForSizeOf_SgTypeDouble(&ad_tyc_1,1);
    FW_ADMM_Allocate(&ad_var_1,ad_tyc_1,0);
    allocateBytesForSizeOf_SgTypeDouble(&ad_tyc_0,1);
    FW_ADMM_Allocate(&ad_var_0,ad_tyc_0,0);
    ad_Symbol_4 = 0;
    for (i = 1; i <= n - 2; i = i + 1) {
      if (0 < DERIV_val(u_1[i])) {
        DERIV_val(ad_var_0) = 0;
        ad_Symbol_5 = 1;
        push_i_s0(ad_Symbol_5);
      }
       else {
        DERIV_val(ad_var_0) = DERIV_val(u_1[i]);
        push_i_s0(i);
        ad_Symbol_6 = 0;
        push_i_s0(ad_Symbol_6);
      }
      if (0 > DERIV_val(u_1[i])) {
        DERIV_val(ad_var_1) = 0;
        ad_Symbol_7 = 1;
        push_i_s0(ad_Symbol_7);
      }
       else {
        DERIV_val(ad_var_1) = DERIV_val(u_1[i]);
        push_i_s0(i);
        ad_Symbol_8 = 0;
        push_i_s0(ad_Symbol_8);
      }
      push_p_s0(ad_var_2);
      ad_var_2 = &u[i];
      ad_aux_2 = -DERIV_val( *C);
      ad_aux_0 = -DERIV_val(u_1[i]) + DERIV_val(u_1[i + 1]);
      ad_aux_1 = DERIV_val(u_1[i]) - DERIV_val(u_1[i - 1]);
      ad_aux_3 = ad_aux_0 * DERIV_val(ad_var_0) + ad_aux_1 * DERIV_val(ad_var_1);
      ad_aux_4 = - 2.00000;
      ad_aux_5 = ad_aux_4 * DERIV_val(u_1[i]) + DERIV_val(u_1[i - 1]) + DERIV_val(u_1[i + 1]);
      ad_lin_4 = ad_aux_3;
      ad_lin_0 = DERIV_val(ad_var_0);
      ad_lin_1 = ad_aux_0;
      ad_lin_2 = DERIV_val(ad_var_1);
      ad_lin_3 = ad_aux_1;
      ad_lin_5 = ad_aux_2;
      ad_lin_7 = ad_aux_5;
      ad_lin_6 = ad_aux_4;
      ad_lin_8 = DERIV_val( *D);
      DERIV_val(ad_var_3) = DERIV_val( *ad_var_2) + (ad_aux_2 * ad_aux_3 + DERIV_val( *D) * ad_aux_5 + DERIV_val(u_1[i]));
      ad_acc_0 = -1 * ad_lin_4;
      ad_Symbol_18 = i + 1;
      push_i_s0(ad_Symbol_18);
      ad_Symbol_19 = i - 1;
      push_i_s0(ad_Symbol_19);
      ad_Symbol_20 = i - 1;
      push_i_s0(ad_Symbol_20);
      ad_Symbol_21 = i + 1;
      push_i_s0(ad_Symbol_21);
      push_s0(ad_lin_0);
      push_s0(ad_lin_1);
      push_s0(ad_lin_2);
      push_s0(ad_lin_3);
      push_s0(ad_lin_6);
      push_s0(ad_lin_5);
      push_s0(ad_acc_0);
      push_s0(ad_lin_7);
      push_s0(ad_lin_8);
      push_i_s0(i);
      DERIV_val( *ad_var_2) = DERIV_val(ad_var_3);
      ad_Symbol_4 = ad_Symbol_4 + 1;
    }
    push_i_s0(ad_Symbol_4);
    push_p_s0(ad_var_2);
    FW_ADMM_Deallocate(&ad_var_0,0,0);
    FW_ADMM_Deallocate(&ad_var_1,0,0);
    FW_ADMM_Deallocate(&ad_var_3,0,0);
    FW_ADMM_Deallocate(&i,0,0);
  }
   else if (our_rev_mode . adjoint == 1) {
    BW_ADMM_Deallocate(&i,0,0,0);
    BW_ADMM_Deallocate(&ad_var_3,0,0,0);
    BW_ADMM_Deallocate(&ad_var_1,0,0,0);
    BW_ADMM_Deallocate(&ad_var_0,0,0,0);
    pop_p_s0(ad_var_2);
    ADMM_Rebase(&ad_var_2);
    pop_i_s0(ad_Symbol_0);
    for (ad_Symbol_1 = 1; ad_Symbol_1 <= ad_Symbol_0; ad_Symbol_1 = ad_Symbol_1 + 1) {
      IncDeriv(ad_prp_17, *ad_var_2);
      ZeroDeriv( *ad_var_2);
      IncDeriv(ad_var_3,ad_prp_17);
      ZeroDeriv(ad_prp_17);
      pop_i_s0(ad_Symbol_22);
      pop_s0(ad_Symbol_23);
      pop_s0(ad_Symbol_24);
      pop_s0(ad_Symbol_25);
      pop_s0(ad_Symbol_26);
      pop_s0(ad_Symbol_27);
      pop_s0(ad_Symbol_28);
      pop_s0(ad_Symbol_29);
      pop_s0(ad_Symbol_30);
      pop_s0(ad_Symbol_31);
      Saxpy(ad_Symbol_23,ad_var_3,ad_prp_16);
      Saxpy(ad_Symbol_24,ad_var_3,ad_prp_8);
      Saxpy(ad_Symbol_25,ad_var_3,ad_prp_1);
      Saxpy(ad_Symbol_26,ad_var_3,ad_prp_15);
      IncDeriv(ad_prp_12,ad_var_3);
      IncDeriv(ad_prp_0,ad_var_3);
      ZeroDeriv(ad_var_3);
      Saxpy(ad_Symbol_27,ad_prp_16,ad_prp_9);
      IncDeriv(ad_prp_10,ad_prp_16);
      IncDeriv(ad_prp_11,ad_prp_16);
      ZeroDeriv(ad_prp_16);
      Saxpy(ad_Symbol_28,ad_prp_15,ad_prp_7);
      Saxpy(ad_Symbol_29,ad_prp_15,ad_prp_14);
      Saxpy(ad_Symbol_30,ad_prp_15,ad_prp_4);
      Saxpy(ad_Symbol_31,ad_prp_15,ad_prp_13);
      ZeroDeriv(ad_prp_15);
      DecDeriv(ad_prp_6,ad_prp_14);
      IncDeriv(ad_prp_5,ad_prp_14);
      ZeroDeriv(ad_prp_14);
      DecDeriv(ad_prp_2,ad_prp_13);
      IncDeriv(ad_prp_3,ad_prp_13);
      ZeroDeriv(ad_prp_13);
      IncDeriv(u_1[ad_Symbol_22],ad_prp_12);
      ZeroDeriv(ad_prp_12);
      pop_i_s0(ad_Symbol_32);
      IncDeriv(u_1[ad_Symbol_32],ad_prp_11);
      ZeroDeriv(ad_prp_11);
      pop_i_s0(ad_Symbol_33);
      IncDeriv(u_1[ad_Symbol_33],ad_prp_10);
      ZeroDeriv(ad_prp_10);
      IncDeriv(u_1[ad_Symbol_22],ad_prp_9);
      ZeroDeriv(ad_prp_9);
      IncDeriv( *D,ad_prp_8);
      ZeroDeriv(ad_prp_8);
      IncDeriv(ad_var_1,ad_prp_7);
      ZeroDeriv(ad_prp_7);
      pop_i_s0(ad_Symbol_34);
      IncDeriv(u_1[ad_Symbol_34],ad_prp_6);
      ZeroDeriv(ad_prp_6);
      IncDeriv(u_1[ad_Symbol_22],ad_prp_5);
      ZeroDeriv(ad_prp_5);
      IncDeriv(ad_var_0,ad_prp_4);
      ZeroDeriv(ad_prp_4);
      pop_i_s0(ad_Symbol_35);
      IncDeriv(u_1[ad_Symbol_35],ad_prp_3);
      ZeroDeriv(ad_prp_3);
      IncDeriv(u_1[ad_Symbol_22],ad_prp_2);
      ZeroDeriv(ad_prp_2);
      IncDeriv( *C,ad_prp_1);
      ZeroDeriv(ad_prp_1);
      IncDeriv( *ad_var_2,ad_prp_0);
      ZeroDeriv(ad_prp_0);
      pop_p_s0(ad_var_2);
      ADMM_Rebase(&ad_var_2);
      pop_i_s0(ad_Symbol_2);
      if (ad_Symbol_2 != 0) {
        ZeroDeriv(ad_var_1);
      }
       else {
        pop_i_s0(ad_Symbol_36);
        IncDeriv(ad_prp_18,ad_var_1);
        ZeroDeriv(ad_var_1);
        IncDeriv(u_1[ad_Symbol_36],ad_prp_18);
        ZeroDeriv(ad_prp_18);
      }
      pop_i_s0(ad_Symbol_3);
      if (ad_Symbol_3 != 0) {
        ZeroDeriv(ad_var_0);
      }
       else {
        pop_i_s0(ad_Symbol_37);
        IncDeriv(ad_prp_19,ad_var_0);
        ZeroDeriv(ad_var_0);
        IncDeriv(u_1[ad_Symbol_37],ad_prp_19);
        ZeroDeriv(ad_prp_19);
      }
    }
    BW_ADMM_Allocate(&ad_var_0,0);
    BW_ADMM_Allocate(&ad_var_1,0);
    BW_ADMM_Allocate(&ad_var_3,0);
    BW_ADMM_Allocate(&i,0);
  }
}
