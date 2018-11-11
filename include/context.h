/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file control/context.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2018-11-11
 *
 **/
#ifndef _MORSE_CONTEXT_H_
#define _MORSE_CONTEXT_H_

#include "chameleon/morse_struct.h"

/*******************************************************************************
 *  Routines to handle threads context
 **/
#ifdef __cplusplus
extern "C" {
#endif

MORSE_context_t* morse_context_create  ();
MORSE_context_t* morse_context_self    ();
int              morse_context_destroy ();

#ifdef __cplusplus
}
#endif

#endif
