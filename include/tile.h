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
 * @file tile.h
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
#ifndef _MORSE_TILE_H_
#define _MORSE_TILE_H_

#if defined( _WIN32 ) || defined( _WIN64 )
typedef __int64 int64_t;
#else
#include <inttypes.h>
#endif

#define ELTADDR(A, type, m, n)  (type *)morse_geteltaddr(A, m, n)
#define ELTLDD(A, k) ( ( (((k)-1)/(A).mb) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *  Internal routines - dynamic scheduling
 **/
void morse_pztile_to_lapack(MORSE_desc_t*, MORSE_Complex64_t*, int, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctile_to_lapack(MORSE_desc_t*, MORSE_Complex32_t*, int, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtile_to_lapack(MORSE_desc_t*, double*, int, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstile_to_lapack(MORSE_desc_t*, float*, int, MORSE_sequence_t *sequence, MORSE_request_t *request);

void morse_pzlapack_to_tile(MORSE_Complex64_t*, int, MORSE_desc_t*, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclapack_to_tile(MORSE_Complex32_t*, int, MORSE_desc_t*, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlapack_to_tile(double*, int, MORSE_desc_t*, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslapack_to_tile(float*, int, MORSE_desc_t*, MORSE_sequence_t *sequence, MORSE_request_t *request);

void morse_pztile_zero(MORSE_desc_t *dA, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctile_zero(MORSE_desc_t *dA, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtile_zero(MORSE_desc_t *dA, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstile_zero(MORSE_desc_t *dA, MORSE_sequence_t *sequence, MORSE_request_t *request);

#ifdef __cplusplus
}
#endif

#endif
