/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dgemv.c
 *
 * StarPU codelet to compute dense matrix-vector multiplication.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include <cblas.h>
#include "../include/starpu_exageostat.h"

static void CORE_dgemv_starpu(void *buffers[], void *cl_arg) {
    int m, n, ldam;
    int m0, n0;
    double* A;
    double* Z;
    double* Zout;

    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    Z = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    Zout = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &ldam, &m0, &n0);
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1, A, m, Z, 1, 0, Zout, 1);
}

static struct starpu_codelet cl_dgemv =
        {
                .where         = STARPU_CPU,
                .cpu_funcs     = {CORE_dgemv_starpu},
                .nbuffers     = 3,
                .modes         = {STARPU_R, STARPU_R, STARPU_RW},
                .name         = "dgemv"
        };

//*******************************************************************************
static void CORE_sgemv_starpu(void *buffers[], void *cl_arg) {
    int m, n, ldam;
    int m0, n0;
    float *A;
    float *Z;
    float *Zout;

    A = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);
    Z = (float *) STARPU_MATRIX_GET_PTR(buffers[1]);
    Zout = (float *) STARPU_MATRIX_GET_PTR(buffers[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &ldam, &m0, &n0);
    cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, 1, A, m, Z, 1, 0, Zout, 1);
}

static struct starpu_codelet cl_sgemv =
        {
                .where           = STARPU_CPU,
                .cpu_funcs       = {CORE_sgemv_starpu},
                .nbuffers        = 3,
                .modes           = {STARPU_R, STARPU_R, STARPU_RW},
                .name            = "sgemv"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_dgemv_Tile_Async - Matrix Vector Multiplication.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           Chameleon descriptor.
 *
 * @param[in] descZ
 *           Chameleon descriptor.
 *
 * @param[out] descZout
 *           Chameleon descriptor.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[in] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int EXAGEOSTAT_MLE_dgemv_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descZ, CHAM_desc_t *descZout,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, n, ldam, m0, n0;
    int tempmm, tempnn;
    CHAM_desc_t A = *descA;

    struct starpu_codelet *cl = &cl_dgemv;
    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            ldam = mBLKLDD(descA, m);
            m0 = m * A.mb;
            n0 = n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &ldam, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, n),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descZ, ChamRealDouble, m, 0),
                               STARPU_RW, EXAGEOSTAT_RTBLKADDR(descZout, ChamRealDouble, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                               STARPU_NAME, "dgemv",
#endif
                               0);
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    CHAMELEON_Sequence_Wait(sequence);
    return CHAMELEON_SUCCESS;
}

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex32_t_Tile (Single precision).
 *
 *  EXAGEOSTAT_MLE_sgemv_Tile_Async - Matrix Vector Multiplication.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           Chameleon descriptor.
 *
 * @param[in] descZ
 *           Chameleon descriptor.
 *
 * @param[out] descZout
 *           Chameleon descriptor.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[in] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int EXAGEOSTAT_MLE_sgemv_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descZ, CHAM_desc_t *descZout,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, n, ldam, m0, n0;
    int tempmm, tempnn;
    CHAM_desc_t A = *descA;

    struct starpu_codelet *cl = &cl_sgemv;
    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            ldam = mBLKLDD(descA, m);
            m0 = m * A.mb;
            n0 = n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &ldam, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealFloat, A.mt, n),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descZ, ChamRealFloat, 0, n),
                               STARPU_RW, EXAGEOSTAT_RTBLKADDR(descZout, ChamRealFloat, 0, n),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                               STARPU_NAME, "sgemv",
#endif
                               0);
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    CHAMELEON_Sequence_Wait(sequence);
    return CHAMELEON_SUCCESS;
}