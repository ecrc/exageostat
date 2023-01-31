/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_ddotp.c
 *
 * StarPU codelet to Calculate the dot product of the Z vector.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_ddotp_starpu(void *buffers[], void *cl_arg) {
    int m, m0;
    int transA, transB;
    int index;
    double* A;
    double* B;
    double* C;
    int increment;

    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    B = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    C = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &m0, &transA, &transB, &index, &increment);

    int i = 0;
    double local_dot = 0;
    for (i = index; i < m; i += increment) {
        local_dot += A[i] * B[i];
    }
    C[0] += local_dot;
}

static struct starpu_codelet cl_ddotp =
        {
                .where          = STARPU_CPU,
                .cpu_funcs      = {CORE_ddotp_starpu},
                .nbuffers       = 3,
                .modes          = {STARPU_R, STARPU_R, STARPU_RW},
                .name           = "ddotp"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_ddotp_Async- codelet to compute dot product of A.A.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           Chameleon descriptor
 *
 * @param[out] descproduct
 *           dot product descriptor.
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
int EXAGEOSTAT_MLE_ddotp_Async(CHAM_enum transA, CHAM_enum transB,
                              CHAM_desc_t *descA, CHAM_desc_t *descB,
                              CHAM_desc_t *descC, CHAM_enum index, CHAM_enum increment,
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
    struct starpu_codelet *cl = &cl_ddotp;

    printf("%d - %d -%d\n", descA->m, descB->m, descC->m);
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &transA, sizeof(CHAM_enum),
                           STARPU_VALUE, &transB, sizeof(CHAM_enum),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, 0),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descB, ChamRealDouble, m, 0),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descC, ChamRealDouble, 0, 0),
                           STARPU_VALUE, &index, sizeof(CHAM_enum),
                           STARPU_VALUE, &increment, sizeof(CHAM_enum),
                           0);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}