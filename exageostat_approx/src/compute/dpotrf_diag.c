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
 *  CHAMELEON computational routines
 *  CHAMELEON is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2022-11-09
 * @generated d Fri Dec  1 14:38:30 2017
 *
 **/
#include "../include/diag.h"
#include <chameleon/constants.h>
#include <global.h>
#include <common.h>
#include <compute_d.h>

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  CHAMELEON_dpotrf - Computes the Cholesky factorization of a symmetric positive definite
 *  (or Symmetric positive definite in the real case) matrix A.
 *  The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = ChamLower}^{U^H\times U, if uplo = ChamUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Symmetric) matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
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
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_dpotrf_Tile
 * @sa CHAMELEON_dpotrf_Tile_Async
 * @sa CHAM_cpotrf
 * @sa CHAM_dpotrf
 * @sa CHAM_spotrf
 * @sa CHAM_dpotrs
 *
 ******************************************************************************/
int CHAM_dpotrf_diag(CHAM_enum uplo, int N,
                     double* A, int LDA, int diag_thick) {
    int NB;
    int status;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t *request = CHAMELEON_SUCCESS;
    CHAM_desc_t descAl, descAt;

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAM_diag_dpotrf", "CHAM not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        chameleon_error("CHAM_diag_dpotrf", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAM_diag_dpotrf", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, N)) {
        chameleon_error("CHAM_diag_dpotrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_max(N, 0) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_DPOSV, N, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAM_diag_dpotrf", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create(chamctxt, &sequence);

    /* Submit the matrix conversion */
    chameleon_dlap2tile(chamctxt, &descAl, &descAt, ChamDescInout, uplo,
                        A, NB, NB, LDA, N, N, N, sequence, &request);

    /* Call the tile interface */
    CHAM_dpotrf_diag_Tile_Async(uplo, &descAt, diag_thick, sequence, &request);

    /* Submit the matrix conversion back */
    chameleon_dtile2lap(chamctxt, &descAl, &descAt,
                        ChamDescInout, uplo, sequence, &request);

    chameleon_sequence_wait(chamctxt, sequence);

    /* Cleanup the temporary data */
    chameleon_dtile2lap_cleanup(chamctxt, &descAl, &descAt);

    status = sequence->status;
    chameleon_sequence_destroy(chamctxt, sequence);

    return status;


}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  CHAMELEON_dpotrf_Tile - Computes the Cholesky factorization of a symmetric positive definite
 *  or Symmetric positive definite matrix.
 *  Tile equivalent of CHAM_dpotrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A is stored;
 *          = ChamLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite (or Symmetric) matrix A.
 *          If uplo = ChamUpper, the leading N-by-N upper triangular part of A
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
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa CHAM_dpotrf
 * @sa CHAMELEON_dpotrf_Tile_Async
 * @sa CHAM_cpotrf_Tile
 * @sa CHAMELEON_dpotrf_Tile
 * @sa CHAMELEON_spotrf_Tile
 * @sa CHAM_dpotrs_Tile
 *
 ******************************************************************************/
int CHAM_dpotrf_diag_Tile(CHAM_enum uplo, CHAM_desc_t *A, int diag_thick) {
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t *request = CHAMELEON_SUCCESS;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile", "CHAM not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create(chamctxt, &sequence);
    CHAM_dpotrf_diag_Tile_Async(uplo, A, diag_thick, sequence, &request);
    CHAMELEON_Desc_Flush(A, sequence);
    chameleon_sequence_wait(chamctxt, sequence);
    status = sequence->status;
    chameleon_sequence_destroy(chamctxt, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  CHAMELEON_dpotrf_Tile_Async - Computes the Cholesky factorization of a symmetric
 *  positive definite or Symmetric positive definite matrix.
 *  Non-blocking equivalent of CHAMELEON_dpotrf_Tile().
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
 * @sa CHAM_dpotrf
 * @sa CHAMELEON_dpotrf_Tile
 * @sa CHAM_cpotrf_Tile_Async
 * @sa CHAMELEON_dpotrf_Tile_Async
 * @sa CHAMELEON_spotrf_Tile_Async
 * @sa CHAM_dpotrs_Tile_Async
 *
 ******************************************************************************/
int CHAM_dpotrf_diag_Tile_Async(CHAM_enum uplo, CHAM_desc_t *A, int diag_thick,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile_Async", "CHAM not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_dpotrf_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS)
        request->status = CHAMELEON_SUCCESS;
    else
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_dpotrf_Tile_Async", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        chameleon_error("CHAMELEON_dpotrf_Tile_Async", "only square tiles supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    if (uplo != ChamUpper && uplo != ChamLower) {
        chameleon_error("CHAMELEON_dpotrf_Tile_Async", "illegal value of uplo");
        return chameleon_request_fail(sequence, request, -1);
    }
    /* Quick return */
    CHAM_pdpotrf_diag(uplo, A, diag_thick, sequence, request);

    return CHAMELEON_SUCCESS;
}
