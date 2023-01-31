/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dmdet.c
 *
 * StarPU codelet to Calculate determinant of a given triangular matrix (A)
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_dmdet_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    double* A;
    int m0;
    int n0;
    double det = 0;
    double* determinant = &det;

    *determinant = 0;
    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    determinant = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);
    double local_det = core_dmdet(A, m, n, m0, n0);
    *determinant += local_det;
}

static struct starpu_codelet cl_dmdet =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_dmdet_starpu},
                .nbuffers    = 2,
                .modes        = {STARPU_R, STARPU_RW},
                .name        = "dmdet"
        };

static void CORE_smdet_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    float *A;
    int m0;
    int n0;
    double det = 0;
    double* determinant = &det;

    *determinant = 0;
    A = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);
    determinant = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);
    float local_det = core_smdet(A, m, n, m0, n0);
    *determinant += local_det;
}

static struct starpu_codelet cl_smdet =
        {
                .where          = STARPU_CPU,
                .cpu_funcs      = {CORE_smdet_starpu},
                .nbuffers       = 2,
                .modes          = {STARPU_R, STARPU_RW},
                .name           = "smdet"
        };

/***************************************************************************//**
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_dmdet_Tile_Async  - Calculate determinant for triangular matrix.
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
 * @param[out] descdet
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
int EXAGEOSTAT_MLE_dmdet_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request,
                                   CHAM_desc_t *descdet) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0, n0;
    int tempmm;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_dmdet;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, m),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descdet, ChamRealDouble, 0, 0),
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
 *  EXAGEOSTAT_MLE_smdet_Tile_Async  - Calculate determinant for triangular matrix.
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
 * @param[out] descdet
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
int EXAGEOSTAT_MLE_smdet_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request,
                                   CHAM_desc_t *descdet) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0, n0;
    int tempmm;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_smdet;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealFloat, m, m),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descdet, ChamRealFloat, 0, 0),
                           0);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}