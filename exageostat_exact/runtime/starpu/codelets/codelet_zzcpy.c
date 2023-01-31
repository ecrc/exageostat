/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dzcpy.c
 *
 * StarPU codelet to Copy contents of descriptor to vector
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include "../include/starpu_exageostat.h"
#include "exageostatcore.h"

static void CORE_dzcpy_starpu(void *buffers[], void *cl_arg) {
    int m;
    double* A;
    int m0;
    double* r;

    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &m0, &r);
    core_dzcpy(A, m, m0, r);
}

static struct starpu_codelet cl_dzcpy =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_dzcpy_starpu},
                .nbuffers    = 1,
                .modes        = {STARPU_W},
                .name        = "dzcpy"
        };


static void CORE_szcpy_starpu(void *buffers[], void *cl_arg) {
    int m;
    float *A;
    int m0;
    float *r;

    A = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &m0, &r);

    core_szcpy(A, m, m0, r);
}

static struct starpu_codelet cl_szcpy =
        {
                .where          = STARPU_CPU,
                .cpu_funcs      = {CORE_szcpy_starpu},
                .nbuffers       = 1,
                .modes          = {STARPU_W},
                .name           = "szcpy"
        };

/***************************************************************************//**
 *
 * @ingroup CHAMELEON_Complex64_t_Tile (double precision).
 *
 *  CHAMELEON_MLE_zcpy_Tile_Async - copy Chameleon descriptor to vector dobule *.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[out] descA
 *           Chameleon descriptor
 *
 * @param[in] r
 *           double*
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
int EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAM_desc_t *descA, double* r, RUNTIME_sequence_t *sequence,
                                   RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_dzcpy;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;

        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &r, sizeof(double),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                           STARPU_NAME, "dzcpy",
#endif
                           0);

    }
    RUNTIME_options_ws_free(&options);
    return CHAMELEON_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup CHAMELEON_Complex32_t_Tile (single precision).
 *
 *  EXAGEOSTAT_MLE_szcpy_Tile_Async - copy Chameleon descriptor to vector float*.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[out] descA
 *           Chameleon descriptor
 *
 * @param[in] r
 *           double*
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
int EXAGEOSTAT_MLE_szcpy_Tile_Async(CHAM_desc_t *descA, double* r, RUNTIME_sequence_t *sequence,
                                   RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_szcpy;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;

        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &r, sizeof(double),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealFloat, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                           STARPU_NAME, "szcpy",
#endif
                           0);

    }
    RUNTIME_options_ws_free(&options);
    return CHAMELEON_SUCCESS;
}