
#include "ad_rev.h"
     
__ADIC_RevType our_rev_mode = {0,0,0,0,0,0,0};
//__ADIC_RevType our_rev_mode = {0};

/*
int isRevModePlain(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.plain;
}

int isRevModeTape(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.tape;
}

int isRevModeAdjoint(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.adjoint;
}

int isRevModeArgStore(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.arg_store;
}

int isRevModeArgRestore(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.arg_restore;
}

int isRevModeRedStore(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.res_store;
}

int isRevModeResRestore(__ADIC_RevType our_rev_mode)
{
  return  our_rev_mode.res_restore;
}

__ADIC_RevType getRevMode()
{
  return our_rev_mode;  
}*/

void __ADIC_RevInit(__ADIC_RevType *adic_rev)
     {
      adic_rev->arg_store = FALSE;
      adic_rev->arg_restore = FALSE;
      adic_rev->res_store = FALSE;
      adic_rev->res_restore = FALSE;
      adic_rev->plain = TRUE;
      adic_rev->tape = FALSE;
      adic_rev->adjoint = FALSE;
     }
     
     void __ADIC_RevPlain(__ADIC_RevType *adic_rev)
     {
      adic_rev->arg_store = FALSE;
      adic_rev->arg_restore = FALSE;
      adic_rev->res_store = FALSE;
      adic_rev->res_restore = FALSE;
      adic_rev->plain = TRUE;
      adic_rev->tape = FALSE;
      adic_rev->adjoint = FALSE;
     }
     void __ADIC_RevTape(__ADIC_RevType *adic_rev)
     {
      adic_rev->arg_store = FALSE;
      adic_rev->arg_restore = FALSE;
      adic_rev->res_store = FALSE;
      adic_rev->res_restore = FALSE;
      adic_rev->plain = FALSE;
      adic_rev->tape = TRUE;
      adic_rev->adjoint = FALSE;
     }
     void __ADIC_RevAdjoint(__ADIC_RevType *adic_rev)
     {
      adic_rev->arg_store = FALSE;
      adic_rev->arg_restore = FALSE;
      adic_rev->res_store = FALSE;
      adic_rev->res_restore = FALSE;
      adic_rev->plain = FALSE;
      adic_rev->tape = FALSE;
      adic_rev->adjoint = TRUE;
     }
     void __ADIC_RevStorePlain(__ADIC_RevType *adic_rev)
     {
      adic_rev->arg_store = TRUE;
      adic_rev->arg_restore = FALSE;
      adic_rev->res_store = FALSE;
      adic_rev->res_restore = FALSE;
      adic_rev->plain = TRUE;
      adic_rev->tape = FALSE;
      adic_rev->adjoint = FALSE;
     }
     void __ADIC_RevRestoreTape(__ADIC_RevType *adic_rev)
     {
      adic_rev->arg_store = FALSE;
      adic_rev->arg_restore = TRUE;
      adic_rev->res_store = FALSE;
      adic_rev->res_restore = FALSE;
      adic_rev->plain = FALSE;
      adic_rev->tape = TRUE;
      adic_rev->adjoint = FALSE;
     }

