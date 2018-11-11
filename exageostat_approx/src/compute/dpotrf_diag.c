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
 * @file dpotrf.c
 *
 *  MORSE computational routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-11
 * @generated d Fri Dec  1 14:38:30 2017
 *
 **/
#include "../include/diag.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  MORSE_dpotrf - Computes the Cholesky factorization of a symmetric positive definite
 *  (or Symmetric positive definite in the real case) matrix A.
 *  The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = MorseLower}^{U^H\times U, if uplo = MorseUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Symmetric) matrix A.
 *          If uplo = MorseUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa MORSE_dpotrf_Tile
 * @sa MORSE_dpotrf_Tile_Async
 * @sa MORSE_cpotrf
 * @sa MORSE_dpotrf
 * @sa MORSE_spotrf
 * @sa MORSE_dpotrs
 *
 ******************************************************************************/
int MORSE_dpotrf_diag(MORSE_enum uplo, int N,
                  double *A, int LDA, int diag_thick)
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descAl, descAt;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_diag_dpotrf", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_diag_dpotrf", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_diag_dpotrf", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        morse_error("MORSE_diag_dpotrf", "illegal value of LDA");
        return -4;
    }
/* Quick return */
    if (chameleon_max(N, 0) == 0)
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_DPOSV, N, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_diag_dpotrf", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_dlap2tile( morse, &descAl, &descAt, MorseDescInout, uplo,
                     A, NB, NB, LDA, N, N, N, sequence, &request );

    /* Call the tile interface */
    MORSE_dpotrf_diag_Tile_Async(uplo, &descAt, diag_thick, sequence, &request);

    /* Submit the matrix conversion back */
    morse_dtile2lap( morse, &descAl, &descAt,
                     MorseDescInout, uplo, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_dtile2lap_cleanup( morse, &descAl, &descAt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );

    return status;


}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  MORSE_dpotrf_Tile - Computes the Cholesky factorization of a symmetric positive definite
 *  or Symmetric positive definite matrix.
 *  Tile equivalent of MORSE_dpotrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite (or Symmetric) matrix A.
 *          If uplo = MorseUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa MORSE_dpotrf
 * @sa MORSE_dpotrf_Tile_Async
 * @sa MORSE_cpotrf_Tile
 * @sa MORSE_dpotrf_Tile
 * @sa MORSE_spotrf_Tile
 * @sa MORSE_dpotrs_Tile
 *
 ******************************************************************************/
int MORSE_dpotrf_diag_Tile(MORSE_enum uplo, MORSE_desc_t *A, int diag_thick)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_dpotrf_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    MORSE_dpotrf_diag_Tile_Async(uplo, A, diag_thick, sequence, &request);
    MORSE_Desc_Flush( A, sequence );
    morse_sequence_wait(morse, sequence);
   // RUNTIME_desc_getoncpu(A);

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  MORSE_dpotrf_Tile_Async - Computes the Cholesky factorization of a symmetric
 *  positive definite or Symmetric positive definite matrix.
 *  Non-blocking equivalent of MORSE_dpotrf_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa MORSE_dpotrf
 * @sa MORSE_dpotrf_Tile
 * @sa MORSE_cpotrf_Tile_Async
 * @sa MORSE_dpotrf_Tile_Async
 * @sa MORSE_spotrf_Tile_Async
 * @sa MORSE_dpotrs_Tile_Async
 *
 ******************************************************************************/
int MORSE_dpotrf_diag_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, int diag_thick,
                             MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_dpotrf_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_dpotrf_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_dpotrf_Tile_Async", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_dpotrf_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_dpotrf_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (uplo != MorseUpper && uplo != MorseLower) {
        morse_error("MORSE_dpotrf_Tile_Async", "illegal value of uplo");
        return morse_request_fail(sequence, request, -1);
    }
    /* Quick return */
/*
    if (chameleon_max(N, 0) == 0)
        return MORSE_SUCCESS;
*/
    morse_pdpotrf_diag(uplo, A, diag_thick, sequence, request);

    return MORSE_SUCCESS;
}
