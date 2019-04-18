/*
 * ad_grad.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: snarayan
 */

#ifndef AD_GRAD_HPP_
#define AD_GRAD_HPP_
#include <stdio.h>
#include <stdlib.h>
#include "ad_types.h"

void
ADIC_Init();

void
ADIC_Finalize();

void
ADIC_SetForwardMode();

void
ADIC_SetReverseMode();

int
ADIC_GetMode();

/**
 * Internal ADIC macros and functions, do not use in directly.
 */
int
__ADIC_IncrShadowVar();
void
__ADIC_CommitShadowVar();
void
__ADIC_ResetShadowVar();
int
__ADIC_GradSize();

extern int __adic_grad_size;
#define __ADIC_GradSize() (__adic_grad_size)

#define __ADIC_Validate(g)

/* Setting independent and dependent variables
 * Prototype: void ADIC_SetIndep(DERIV_TYPE var);
 */
#define ADIC_SetIndep(var)                                              \
  {                                                                     \
   int mode = ADIC_GetMode();                                          \
   int __adic_pos = -1;\
    if(mode==0) {                                                  \
      __adic_pos = __ADIC_IncrShadowVar();                            \
      if (__adic_pos > ADIC_GRADVEC_LENGTH) {                             \
        fprintf(stderr, "ADIC_SetIndep() Error: the number of independent " \
                "variables %d exceeds the max (ADIC_GRADVEC_LENGTH=%d)!\n", \
                __adic_pos, ADIC_GRADVEC_LENGTH);                         \
        fflush(stderr);                                                   \
        abort();                                                          \
      }\
      __ADIC_SetGradVal(&(var), 1.0, __adic_pos);                         \
    } else if(mode==1) {                                             \
      __ADIC_SetGradVal(&(var), 0.0, __adic_pos);                         \
    } else                                                                \
      fprintf(stderr, "ADIC_SetIndep() Error: Differentiation Mode not set.\n");\
  }

/* Setting independent and dependent 1-D array variables
 * Prototype: void ADIC_SetIndepArray(DERIV_TYPE arrayvar[], arraysize);
 */
#define ADIC_SetIndepArray(arrayvar, arraysize)                         \
  {                                                                     \
    int __adic_iCtr;                                                    \
    for (__adic_iCtr = 0; __adic_iCtr < arraysize; __adic_iCtr++) {     \
      ADIC_SetIndep((arrayvar)[__adic_iCtr]);                           \
    }                                                                   \
  }

/* Setting independent and dependent variables
 * Prototype: void ADIC_SetIndep(DERIV_TYPE var);
 */

#define ADIC_SetDep(var)                                              \
  {                                                                     \
    int __adic_pos = -1;                            \
    int mode = ADIC_GetMode();                                          \
   if(mode==0) {                                                  \
     __ADIC_SetGradVal(&(var), 0.0, __adic_pos);                         \
   } else if(mode==1) {                                             \
     __adic_pos = __ADIC_IncrShadowVar();                            \
     if (__adic_pos > ADIC_GRADVEC_LENGTH) {                             \
       fprintf(stderr, "ADIC_SetIndep() Error: the number of independent " \
               "variables %d exceeds the max (ADIC_GRADVEC_LENGTH=%d)!\n", \
               __adic_pos, ADIC_GRADVEC_LENGTH);                         \
       fflush(stderr);                                                   \
       abort();                                                          \
     }\
     __ADIC_SetGradVal(&(var), 1.0, __adic_pos);                         \
   } else                                                                \
     fprintf(stderr, "ADIC_SetDep() Error: Differentiation Mode not set.\n");\
  }

/* Setting independent and dependent 1-D array variables
 * Prototype: void ADIC_SetIndepArray(DERIV_TYPE arrayvar[], arraysize);
 */
#define ADIC_SetDepArray(arrayvar, arraysize)                         \
  {                                                                     \
    int __adic_iCtr;                                                    \
    for (__adic_iCtr = 0; __adic_iCtr < arraysize; __adic_iCtr++) {     \
      ADIC_SetDep((arrayvar)[__adic_iCtr]);                           \
    }                                                                   \
  }

#define ADIC_SetIndepDone() __ADIC_CommitShadowVar()
#define ADIC_ResetIndep()   __ADIC_ResetShadowVar()



/* Extract the gradient array
 * Prototype: void ADIC_ExtractGrad(double *tgt, DERIV_TYPE src)
 */
#define ADIC_ExtractGrad(tgt, src)                                      \
  {                                                                     \
    int __adic_iCtr;                                                    \
    for (__adic_iCtr = 0; __adic_iCtr < ADIC_GRADVEC_LENGTH; __adic_iCtr++) { \
      tgt[__adic_iCtr] = DERIV_grad(src)[__adic_iCtr];                  \
    }                                                                   \
  }

/* Extract the value of a DERIV_TYPE variable (same as DERIV_val(src))
 * Prototype: void ADIC_ExtractVal(double *tgt, DERIV_TYPE src)
 */
#define ADIC_ExtractVal(tgt, src)                                      \
  {                                                                    \
    tgt = DERIV_val(src)                                               \
  }

/* Set the value of a DERIV_TYPE variable (same as DERIV_val(src))
 * Prototype: void ADIC_SetValue(DERIV_TYPE tgt, double val)
 */
#define ADIC_SetValue(tgt, val)                                      \
  {                                                                    \
    DERIV_val(tgt) = val;                                               \
  }

/* Set the value of one element of the gradient vector and set the
 * rest of the vector elements to zero.
 * Prototype: void ADIC_SetGradVal(DERIV_TYPE *var, double a, int position);
 */
#define __ADIC_SetGradVal(var, a, i)                                    \
  {                                                                     \
    int __adic_pos2;                                                    \
    __ADIC_Validate(var);                                               \
    for (__adic_pos2 = 0; __adic_pos2 < ADIC_GRADVEC_LENGTH; __adic_pos2++) { \
      if (__adic_pos2 == i) {                                           \
        DERIV_grad(*var)[__adic_pos2] = a;                              \
      } else {                                                            \
        DERIV_grad(*var)[__adic_pos2] = 0.0;                            \
      }                                                                 \
    }                                                                   \
  }

/* Set the value of one element of the gradient vector and set the
 * rest of the vector elements to zero.
 * Prototype: void ADIC_SetGradVal(DERIV_TYPE *var, double a, int position);
 */
#define ADIC_Initialize(var)                                    \
  {                                                        \
  }

#define Saxpy(a, x, y)\
{\
  int __adic_iCtr; double *__adic_gradz = DERIV_grad(y), *__adic_grad1 =  DERIV_grad(x);\
  for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) {\
    __adic_gradz[__adic_iCtr] =  __adic_gradz[__adic_iCtr] + (a) * __adic_grad1[__adic_iCtr];\
  }\
}

#define ADIC_SetDeriv(y,x){\
  int __adic_iCtr; double *__adic_gradz = DERIV_grad(x),\
  *__adic_grads=DERIV_grad(y);\
  for (__adic_iCtr = 0 ; __adic_iCtr < __ADIC_GradSize() ; __adic_iCtr++) {\
    __adic_gradz[__adic_iCtr] = __adic_grads[__adic_iCtr];\
  }\
}

#define DecDeriv(x,y){\
  int __adic_iCtr; double *__adic_gradz = DERIV_grad(x),\
  *__adic_grads=DERIV_grad(y);\
  for (__adic_iCtr = 0 ; __adic_iCtr < __ADIC_GradSize() ; __adic_iCtr++) {\
    __adic_gradz[__adic_iCtr] =  __adic_gradz[__adic_iCtr] - __adic_grads[__adic_iCtr];\
  }\
}

#define IncDeriv(x,y){\
  int __adic_iCtr; double *__adic_gradz = DERIV_grad(x),\
  *__adic_grads=DERIV_grad(y);\
  for (__adic_iCtr = 0 ; __adic_iCtr < __ADIC_GradSize() ; __adic_iCtr++) {\
    __adic_gradz[__adic_iCtr] =  __adic_gradz[__adic_iCtr] + __adic_grads[__adic_iCtr];\
  }\
}

#define ZeroDeriv(x){\
  int __adic_iCtr; double *__adic_gradz = DERIV_grad(x);\
  for (__adic_iCtr = 0 ; __adic_iCtr < __ADIC_GradSize() ; __adic_iCtr++) {\
    __adic_gradz[__adic_iCtr] = 0.0;\
  }\
}

#define ZeroDerivvector(x, __dim0){\
  int __adic_dim0Ctr;\
  for (__adic_dim0Ctr = 0 ;__adic_dim0Ctr < __dim0 ;__adic_dim0Ctr++) {\
    ZeroDeriv((x)[__adic_dim0Ctr])\
  }\
}

#define SetZeroValvector(x, __dim0){\
  int __adic_dim0Ctr;\
  for (__adic_dim0Ctr = 0 ;__adic_dim0Ctr < __dim0 ;__adic_dim0Ctr++) {\
    DERIV_val((x)[__adic_dim0Ctr]) = 0.0;\
  }\
}

#define ZeroDerivmatrix(x, __dim1, __dim0){\
  int __adic_dim0Ctr, __adic_dim1Ctr;\
  for (__adic_dim1Ctr = 0 ;__adic_dim1Ctr < __dim1 ;__adic_dim1Ctr++) {\
    for (__adic_dim0Ctr = 0 ;__adic_dim0Ctr < __dim0 ;__adic_dim0Ctr++) {\
      ZeroDeriv((x)[__adic_dim1Ctr][__adic_dim0Ctr])\
    }\
  }\
}

#define SetZeroValmatrix(x, __dim1, __dim0){\
  int __adic_dim0Ctr, __adic_dim1Ctr;\
  for (__adic_dim1Ctr = 0 ;__adic_dim1Ctr < __dim1 ;__adic_dim1Ctr++) {\
    for (__adic_dim0Ctr = 0 ;__adic_dim0Ctr < __dim0 ;__adic_dim0Ctr++) {\
      DERIV_val((x)[__adic_dim1Ctr][__adic_dim0Ctr]) = 0.0;\
    }\
  }\
}

#define ZeroDerivthree_tensor(x, __dim2, __dim1, __dim0){\
  int __adic_dim0Ctr, __adic_dim1Ctr, __adic_dim2Ctr;\
  for (__adic_dim2Ctr = 0 ;__adic_dim2Ctr < __dim2 ;__adic_dim2Ctr++) {\
    for (__adic_dim1Ctr = 0 ;__adic_dim1Ctr < __dim1 ;__adic_dim1Ctr++) {\
      for (__adic_dim0Ctr = 0 ;__adic_dim0Ctr < __dim0 ;__adic_dim0Ctr++) {\
        ZeroDeriv((x)[__adic_dim2Ctr][__adic_dim1Ctr][__adic_dim0Ctr])\
      }\
    }\
  }\
}

#define SetZeroValthree_tensor(x, __dim2, __dim1, __dim0){\
  int __adic_dim0Ctr, __adic_dim1Ctr, __adic_dim2Ctr;\
  for (__adic_dim2Ctr = 0 ;__adic_dim2Ctr < __dim2 ;__adic_dim2Ctr++) {\
    for (__adic_dim1Ctr = 0 ;__adic_dim1Ctr < __dim1 ;__adic_dim1Ctr++) {\
      for (__adic_dim0Ctr = 0 ;__adic_dim0Ctr < __dim0 ;__adic_dim0Ctr++) {\
        DERIV_val((x)[__adic_dim2Ctr][__adic_dim1Ctr][__adic_dim0Ctr]) = 0.0;\
      }\
    }\
  }\
}

#define MallocDerivPointer(var) {\
  (var) = (DERIV_TYPE *) malloc(sizeof(DERIV_TYPE));\
}

#define oad_AllocateMatching(tgt,src) {\
\
}


#define oad_convert(tgt,src) {\
  (tgt) = (src);\
}

#define oad_ShapeTest(newvar,referencevar) {\
\
}

#endif /* AD_GRAD_HPP_ */
