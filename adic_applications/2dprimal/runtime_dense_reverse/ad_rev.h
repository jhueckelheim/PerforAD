
#ifndef __AD_REV_HPP_
#define __AD_REV_HPP_
  
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

typedef struct modeType
{
     int arg_store;
     int arg_restore;
     int res_store;
     int res_restore;
     int plain;
     int tape;
     int adjoint;
}__ADIC_RevType;

extern __ADIC_RevType our_rev_mode;
  
void __ADIC_RevInit(__ADIC_RevType *adic_rev);

void __ADIC_RevPlain(__ADIC_RevType *adic_rev);

void __ADIC_RevTape(__ADIC_RevType *adic_rev);

void __ADIC_RevAdjoint(__ADIC_RevType *adic_rev);

void __ADIC_RevStorePlain(__ADIC_RevType *adic_rev);

void __ADIC_RevRestoreTape(__ADIC_RevType *adic_rev);
  
/*  int isRevModePlain(__ADIC_RevType our_rev_mode);
  int isRevModeTape(__ADIC_RevType our_rev_mode);
  int isRevModeAdjoint(__ADIC_RevType our_rev_mode);
  int isRevModeArgStore(__ADIC_RevType our_rev_mode);
  int isRevModeArgRestore(__ADIC_RevType our_rev_mode);
  int isRevModeRedStore(__ADIC_RevType our_rev_mode);
  int isRevModeResRestore(__ADIC_RevType our_rev_mode);

  __ADIC_RevType getRevMode();*/
  
#endif

