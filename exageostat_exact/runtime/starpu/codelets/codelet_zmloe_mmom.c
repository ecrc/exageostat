/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dmloe_mmom.c
 *
 * StarPU codelet to Calculate determinant of a given triangular matrix (A)
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-22
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_dmloe_mmom_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    int i;
    double* expr2;
    double* expr3;
    double* expr4;
    double* mloe;
    double* mmom;
    int m0;
    int n0;

    expr2 = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    expr3 = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    expr4 = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
    mloe = (double* ) STARPU_MATRIX_GET_PTR(buffers[3]);
    mmom = (double* ) STARPU_MATRIX_GET_PTR(buffers[4]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);
    double expr2_ = 0, expr3_ = 0, expr4_ = 0;
    printf("--m = %d\n", m);
    for (i = 0; i < m * n; i += 2) {
        expr2_ += expr2[i];
        expr3_ += expr3[i];
        expr4_ += expr4[i];
    }

    if (expr3_ == 0.0){
        *mloe -= 1.0;
    }
    else{
        *mloe += (expr2_ / expr3_) - 1.0;
    }

    if(expr2_ == 0.0){
        *mmom -= 1.0;
    }
    else{
        *mmom += (expr4_ / expr2_) - 1.0;
    }
}

static struct starpu_codelet cl_dmloe_mmom =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_dmloe_mmom_starpu},
                .nbuffers    = 5,
                .modes        = {STARPU_R, STARPU_R, STARPU_R, STARPU_RW, STARPU_RW},
                .name        = "dmloe_mmom"
        };

/***************************************************************************//**
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_dmloe_mmom_Tile_Async  - Calculate determinant for triangular matrix.
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
int EXAGEOSTAT_MLE_dmloe_mmom_Tile_Async(CHAM_desc_t *descexpr2, RUNTIME_sequence_t *descexpr3,
                                        RUNTIME_sequence_t *descexpr4, RUNTIME_sequence_t *descmloe,
                                        RUNTIME_sequence_t *descmmom, RUNTIME_sequence_t *sequence,
                                        RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, n, m0, n0;
    struct starpu_codelet *cl = &cl_dmloe_mmom;
    int tempmm, tempnn;

    for (n = 0; n < descexpr2->nt; n++) {
        tempnn = n == descexpr2->nt - 1 ? descexpr2->n - n * descexpr2->nb : descexpr2->nb;
        for (m = 0; m < descexpr2->mt; m++) {

            tempmm = m == descexpr2->mt - 1 ? descexpr2->m - m * descexpr2->mb : descexpr2->mb;
            m0 = m * descexpr2->mb;
            n0 = n * descexpr2->nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descexpr2, ChamRealDouble, m, n),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descexpr3, ChamRealDouble, m, n),
                               STARPU_R, EXAGEOSTAT_RTBLKADDR(descexpr4, ChamRealDouble, m, n),
                               STARPU_RW, EXAGEOSTAT_RTBLKADDR(descmloe, ChamRealDouble, m, n),
                               STARPU_RW, EXAGEOSTAT_RTBLKADDR(descmmom, ChamRealDouble, m, n),
                               0);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);

    return CHAMELEON_SUCCESS;
}