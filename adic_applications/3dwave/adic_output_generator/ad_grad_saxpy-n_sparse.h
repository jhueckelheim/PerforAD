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
#ifndef __AD_GRAD_SAXPY_N_SPARSE_H
#define __AD_GRAD_SAXPY_N_SPARSE_H

#ifdef ADIC_SPARSE_NO_GRAD
#include "noderiv_sparslinc.h"
#include "noderiv_sparsderiv.hpp"
#endif
#ifdef ADIC_SPARSE
#include "sparslinc.h"
#include "sparsderiv.hpp"
#endif
#ifdef __cplusplus
extern "C" {
#endif
#define ADIC_Sax_1( \
			__adic_a1, __adic_x1, __adic_tgt) {\
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
}

#define ADIC_Sax_2( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, __adic_tgt) {\
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
	std::set<int>::iterator srcSetIter; \
	srcSetIter = (__adic_x2).indexSet.begin(); \
	for (; srcSetIter != (__adic_x2).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
}

#define ADIC_Sax_3( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, __adic_tgt) {\
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
	std::set<int>::iterator srcSetIter; \
	srcSetIter = (__adic_x2).indexSet.begin(); \
	for (; srcSetIter != (__adic_x2).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x3).indexSet.begin(); \
	for (; srcSetIter != (__adic_x3).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
}

#define ADIC_Sax_4( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, __adic_tgt) {\
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
	std::set<int>::iterator srcSetIter; \
	srcSetIter = (__adic_x2).indexSet.begin(); \
	for (; srcSetIter != (__adic_x2).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x3).indexSet.begin(); \
	for (; srcSetIter != (__adic_x3).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x4).indexSet.begin(); \
	for (; srcSetIter != (__adic_x4).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
}

#define ADIC_Sax_5( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, \
			__adic_a5, __adic_x5, __adic_tgt) {\
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
	std::set<int>::iterator srcSetIter; \
	srcSetIter = (__adic_x2).indexSet.begin(); \
	for (; srcSetIter != (__adic_x2).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x3).indexSet.begin(); \
	for (; srcSetIter != (__adic_x3).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x4).indexSet.begin(); \
	for (; srcSetIter != (__adic_x4).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x5).indexSet.begin(); \
	for (; srcSetIter != (__adic_x5).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
}

#define ADIC_Sax_6( \
			__adic_a1, __adic_x1, \
			__adic_a2, __adic_x2, \
			__adic_a3, __adic_x3, \
			__adic_a4, __adic_x4, \
			__adic_a5, __adic_x5, \
			__adic_a6, __adic_x6, __adic_tgt) {\
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
	std::set<int>::iterator srcSetIter; \
	srcSetIter = (__adic_x2).indexSet.begin(); \
	for (; srcSetIter != (__adic_x2).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x3).indexSet.begin(); \
	for (; srcSetIter != (__adic_x3).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x4).indexSet.begin(); \
	for (; srcSetIter != (__adic_x4).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x5).indexSet.begin(); \
	for (; srcSetIter != (__adic_x5).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x6).indexSet.begin(); \
	for (; srcSetIter != (__adic_x6).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
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
	(__adic_tgt).indexSet = (__adic_x1).indexSet; \
	std::set<int>::iterator srcSetIter; \
	srcSetIter = (__adic_x2).indexSet.begin(); \
	for (; srcSetIter != (__adic_x2).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x3).indexSet.begin(); \
	for (; srcSetIter != (__adic_x3).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x4).indexSet.begin(); \
	for (; srcSetIter != (__adic_x4).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x5).indexSet.begin(); \
	for (; srcSetIter != (__adic_x5).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x6).indexSet.begin(); \
	for (; srcSetIter != (__adic_x6).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
	srcSetIter = (__adic_x7).indexSet.begin(); \
	for (; srcSetIter != (__adic_x7).indexSet.end(); srcSetIter++) { \
		(__adic_tgt).indexSet.insert(*srcSetIter); \
	} \
}

#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif
#endif /* #ifndef __AD_GRAD_SAXPY_N_SPARSE_H */
