//#########################################################
// This file is part of OpenAD released under the LGPL.   #
// The full COPYRIGHT notice can be found in the top      #
// level directory of the OpenAD distribution             #
//#########################################################
#include "ad_grad.h"


int __adic_grad_size = 0;
int __adic_total_grad_size = 0;
int __adic_grad_size_shadow = 0;
int __adic_mode = -1;

void ADIC_Init()
{
  __adic_grad_size_shadow = 0;
  __adic_grad_size = 0;
}

void ADIC_Finalize()
{
}

int __ADIC_IncrShadowVar()
{
  return __adic_grad_size_shadow++;
}

void __ADIC_CommitShadowVar()
{
  __adic_grad_size = __adic_grad_size_shadow;
}

void __ADIC_ResetShadowVar()
{
  __adic_grad_size_shadow = 0;
}

void ADIC_SetForwardMode()
{
  __adic_mode = 0;
}

void ADIC_SetReverseMode()
{
  __adic_mode = 1;
}

int ADIC_GetMode()
{
  return __adic_mode;
}
