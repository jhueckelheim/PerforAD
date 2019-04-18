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

/*! \file

 \brief Preaccumulation macros for dense gradient arrays.

 Copyright (c) 2004--2011, UChicago Argonne, LLC <br>
 All rights reserved. <br>
 See $ROSE2XAIF_DIR/Copyright.txt for details. <br>

 */
#ifndef __AD_GRAD_SAXPY_N_DENSE_H
#define __AD_GRAD_SAXPY_N_DENSE_H

#ifndef ADIC_GRAD_LENGTH
#include "ad_types.h"
#else
#include "ad_grad_length_types.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif
#define ADIC_Sax_1( \
			__adic_a1, __adic_x1, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr]; \
	} \
}

#define ADIC_Sax_2( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1), \
			*__adic_grad2 =  DERIV_grad(__adic_x2); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr] + \
			(__adic_a2) * __adic_grad2[__adic_iCtr]; \
	} \
}

#define ADIC_Sax_3( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1), \
			*__adic_grad2 =  DERIV_grad(__adic_x2), \
			*__adic_grad3 =  DERIV_grad(__adic_x3); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr] + \
			(__adic_a2) * __adic_grad2[__adic_iCtr] + \
			(__adic_a3) * __adic_grad3[__adic_iCtr]; \
	} \
}

#define ADIC_Sax_4( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1), \
			*__adic_grad2 =  DERIV_grad(__adic_x2), \
			*__adic_grad3 =  DERIV_grad(__adic_x3), \
			*__adic_grad4 =  DERIV_grad(__adic_x4); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr] + \
			(__adic_a2) * __adic_grad2[__adic_iCtr] + \
			(__adic_a3) * __adic_grad3[__adic_iCtr] + \
			(__adic_a4) * __adic_grad4[__adic_iCtr]; \
	} \
}

#define ADIC_Sax_5( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, \
			__adic_a5, __adic_x5, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1), \
			*__adic_grad2 =  DERIV_grad(__adic_x2), \
			*__adic_grad3 =  DERIV_grad(__adic_x3), \
			*__adic_grad4 =  DERIV_grad(__adic_x4), \
			*__adic_grad5 =  DERIV_grad(__adic_x5); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr] + \
			(__adic_a2) * __adic_grad2[__adic_iCtr] + \
			(__adic_a3) * __adic_grad3[__adic_iCtr] + \
			(__adic_a4) * __adic_grad4[__adic_iCtr] + \
			(__adic_a5) * __adic_grad5[__adic_iCtr]; \
	} \
}

#define ADIC_Sax_6( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, \
			__adic_a5, __adic_x5, \
			__adic_a6, __adic_x6, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1), \
			*__adic_grad2 =  DERIV_grad(__adic_x2), \
			*__adic_grad3 =  DERIV_grad(__adic_x3), \
			*__adic_grad4 =  DERIV_grad(__adic_x4), \
			*__adic_grad5 =  DERIV_grad(__adic_x5), \
			*__adic_grad6 =  DERIV_grad(__adic_x6); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr] + \
			(__adic_a2) * __adic_grad2[__adic_iCtr] + \
			(__adic_a3) * __adic_grad3[__adic_iCtr] + \
			(__adic_a4) * __adic_grad4[__adic_iCtr] + \
			(__adic_a5) * __adic_grad5[__adic_iCtr] + \
			(__adic_a6) * __adic_grad6[__adic_iCtr]; \
	} \
}

#define ADIC_Sax_7( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, \
			__adic_a5, __adic_x5, \
			__adic_a6, __adic_x6, \
			__adic_a7, __adic_x7, __adic_tgt) {\
	int __adic_iCtr; double *__adic_gradz = DERIV_grad(__adic_tgt), \
			*__adic_grad1 =  DERIV_grad(__adic_x1), \
			*__adic_grad2 =  DERIV_grad(__adic_x2), \
			*__adic_grad3 =  DERIV_grad(__adic_x3), \
			*__adic_grad4 =  DERIV_grad(__adic_x4), \
			*__adic_grad5 =  DERIV_grad(__adic_x5), \
			*__adic_grad6 =  DERIV_grad(__adic_x6), \
			*__adic_grad7 =  DERIV_grad(__adic_x7); \
	for (__adic_iCtr = 0; __adic_iCtr < __ADIC_GradSize(); __adic_iCtr++) { \
		__adic_gradz[__adic_iCtr] =  (__adic_a1) * __adic_grad1[__adic_iCtr] + \
			(__adic_a2) * __adic_grad2[__adic_iCtr] + \
			(__adic_a3) * __adic_grad3[__adic_iCtr] + \
			(__adic_a4) * __adic_grad4[__adic_iCtr] + \
			(__adic_a5) * __adic_grad5[__adic_iCtr] + \
			(__adic_a6) * __adic_grad6[__adic_iCtr] + \
			(__adic_a7) * __adic_grad7[__adic_iCtr]; \
	} \
}

#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif
#endif /* #ifndef __AD_GRAD_SAXPY_N_DENSE_H */
