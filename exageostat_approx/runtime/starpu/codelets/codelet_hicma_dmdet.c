/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
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
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"
#include <hicma_struct.h>
#include <hicma_constants.h>
#include "hicma/hicma_ext/control/hicma_context.h"
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
 *           descA:  CHAMELEON descriptor
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 * @param[out] descdet
 *           descerror:  determinant value
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int EXAGEOSTAT_TLR_MLE_dmdet_Tile_Async(HICMA_desc_t *descA, HICMA_sequence_t *sequence, HICMA_request_t *request,
                               HICMA_desc_t *descdet) {

    HICMA_context_t *hicmatxt;
    HICMA_option_t options;
    hicmatxt = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return -2;

    HICMA_RUNTIME_options_init(&options, hicmatxt, sequence, request);
    int m, m0, n0;
    int tempmm;
    HICMA_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_dmdet;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, EXAGEOSTAT_HICMA_RTBLKADDR(descA, HicmaRealDouble, m, 0),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_HICMA_RTBLKADDR(descdet, HicmaRealDouble, 0, 0),
                           0);
    }

    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, hicmatxt);
    return HICMA_SUCCESS;
}
