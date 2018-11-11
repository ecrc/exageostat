/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/***
 *
 *
 * @file async.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2018-11-11
 *
 **/
#ifndef _MORSE_ASYNC_H_
#define _MORSE_ASYNC_H_

#include "chameleon/morse_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *  Internal routines
 **/
int morse_request_fail     (MORSE_sequence_t *sequence, MORSE_request_t *request, int error);
int morse_sequence_create  (MORSE_context_t *MORSE, MORSE_sequence_t **sequence);
int morse_sequence_destroy (MORSE_context_t *MORSE, MORSE_sequence_t *sequence);
int morse_sequence_wait    (MORSE_context_t *MORSE, MORSE_sequence_t *sequence);

#ifdef __cplusplus
}
#endif

#endif
