/*! \file

  \brief This file is for testing all pieces and all steps

  \authors Boyana Norris, Beata Winnicka, Paul Hovland
  \version $Id$

  Copyright (c) 2006, Argonne National Laboratory <br>
  All rights reserved. <br>
  See $ADIC_DIR/share/ADIC_Copyright.txt for details. <br>

*/

#ifndef __AD_TYPES_HPP
#define __AD_TYPES_HPP

#ifdef ADIC_DENSE_REVERSE
#include "adic_gradvec_length.h"
#endif

#ifdef ADIC_DENSE_SEED
#include "adic_gradvec_length_seed.h"
#endif

#include "ad_grad.h"
#include "ad_tape.h"
#include "ad_rev.h"
#include "admm.h"
/*
 * An ultra-simple adic header file for a static-memory vector gradient type.
 */

typedef struct {
	double val;
	double grad[ADIC_GRADVEC_LENGTH];
} DERIV_TYPE;

#define DERIV_val(a)      (a).val
#define DERIV_grad(a)     (a).grad
#define DERIV_gradref(a)  &((a).grad)
#define DERIV_TYPE_ref(a) &(a)

/* Include the first-derivative macros and functions header */

#endif  /* #ifndef __AD_TYPES_HPP */
