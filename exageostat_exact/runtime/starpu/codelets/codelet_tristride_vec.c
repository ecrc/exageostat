/*
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_tristride_vec.c
 *
 * StarPU codelet to Copy contents of descriptor to vector
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2021-01-11
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_tristride_vecstarpu(void *buffers[], void *cl_arg) {
    int m;
    int tempmm;
    double* A;
    double* B;
    double* C;
    double* D;
    int m0;
    int i = 0;
    int j = 0;

    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    B = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    C = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
    D = (double* ) STARPU_MATRIX_GET_PTR(buffers[3]);

    starpu_codelet_unpack_args(cl_arg, &tempmm, &m0, &m);

    //accept only tempmm divided by three (should be optimized)
    j = 0;
    for (i = 0; i < tempmm - 1; i += 3) {
        B[j] = A[i];
        C[j] = A[i + 1];
        D[j] = A[i + 2];
        j++;
    }
}

static struct starpu_codelet cl_tristride_vec =
        {
                .where          = STARPU_CPU,
                .cpu_funcs      = {CORE_tristride_vecstarpu},
                .nbuffers       = 4,
                .modes          = STARPU_R, STARPU_W, STARPU_W, STARPU_W,
                .name           = "tristride_vec"
        };

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
int EXAGEOSTAT_tristride_vec_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descB, CHAM_desc_t *descC, CHAM_desc_t *descD,
                                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t A = *descA;
    CHAM_desc_t B = *descB;
    CHAM_desc_t C = *descC;
    CHAM_desc_t D = *descD;
    struct starpu_codelet *cl = &cl_tristride_vec;
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &m, sizeof(int),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, 0),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(descB, ChamRealDouble, (int) floor(m / 3.0), 0),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(descC, ChamRealDouble, (int) floor(m / 3.0), 0),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(descD, ChamRealDouble, (int) floor(m / 3.0), 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                           STARPU_NAME, "tristride_vec",
#endif
                           0);
    }

    RUNTIME_options_ws_free(&options);
    return CHAMELEON_SUCCESS;
}