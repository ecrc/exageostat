/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dtrace.c
 *
 * StarPU codelet to Calculate trace of a given matrix (A)
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2020-09-18
 *
 **/
#include "../include/starpu_exageostat.h"
#include "exageostatcore.h"

static void CORE_dtrace_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    double* A;
    int m0;
    int n0;
    double s = 0;
    double* sum = &s;
    double* trace;

    *sum = 0;
    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    sum = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    trace = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);

    double local_s = core_dtrace(A, m, n, m0, n0, trace);
    *sum += local_s;
}

static struct starpu_codelet cl_dtrace =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_dtrace_starpu},
                .nbuffers    = 3,
                .modes        = {STARPU_R, STARPU_RW, STARPU_W},
                .name        = "dtrace"
        };


static void CORE_strace_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    float *A;
    int m0;
    int n0;
    double s = 0;
    double* sum = &s;

    *sum = 0;
    A = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);
    sum = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);

    //// TODO: Fix this.
    float local_s = core_strace(A, m, n, m0, n0);
    *sum += local_s;
}

static struct starpu_codelet cl_strace =
        {
                .where          = STARPU_CPU,
                .cpu_funcs      = {CORE_strace_starpu},
                .nbuffers       = 2,
                .modes          = {STARPU_R, STARPU_RW},
                .name           = "strace"
        };

/***************************************************************************//**
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_dtrace_Tile_Async  - Calculate determinant for triangular matrix.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           descA:  Chameleon descriptor
 *
 *
 * @param[out] descsum
 *           descerror:  determinant value
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
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int EXAGEOSTAT_MLE_dtrace_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence,
                                    RUNTIME_request_t *request, CHAM_desc_t *descsum, CHAM_desc_t *desctrace) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0, n0;
    int tempmm;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_dtrace;


    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, m),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descsum, ChamRealDouble, 0, 0),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(desctrace, ChamRealDouble, m, 0),
                           0);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup CHAMELEON_Complex32_t_Tile (single precision)
 *
 *  CHAMELEON_MLE_strace_Tile_Async  - Calculate determinant for triangular matrix.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           descA:  Chameleon descriptor
 *
 *
 * @param[out] descsum
 *           descerror:  determinant value
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
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int CHAMELEON_MLE_strace_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request,
                                    CHAM_desc_t *descsum) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0, n0;
    int tempmm;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_strace;


    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealFloat, m, m),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descsum, ChamRealFloat, 0, 0),
                           0);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}