/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_ng_loglike.c
 *
 * StarPU codelet to Calculate the loglikelihood of non-Gaussian MLE.
 *
 * @version 1.2.0
 *
 * @author Sagnik Mondal
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_ng_loglike_starpu(void *buffers[], void *cl_arg) {
    int m, m0, i;
    double* z;
    double* theta;
    double sum = 0;
    double* s = &sum;

    theta = (double* ) malloc(6 * sizeof(double));
    z = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    s = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);

    starpu_codelet_unpack_args(cl_arg, &m, &m0,
                               &theta[0], &theta[1], &theta[2],
                               &theta[3], &theta[4], &theta[5]);

    double local_sum = core_ng_loglike(z, theta, m);
    *s += local_sum;
}

static struct starpu_codelet cl_ng_loglike =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_ng_loglike_starpu},
                .nbuffers    = 2,
                .modes        = {STARPU_R, STARPU_RW},  //Read access to Z and Read/Write access to the sum.
                .name        = "ng_loglike"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_ng_loglike_Tile_Async - Calculate the loglikelihood of non-Gaussian MLE.
 *  Operates on vectors stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descZ
 *           descZpre:  Observed measurements descZ.
 *
 * @param[out] descSum
 *           descSum:  The loglikelihood Sum of descriptor Z.
 *
 * @param[int]  theta
 *           descerror:  Model parameters.
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

int EXAGEOSTAT_ng_loglike_Tile_Async(CHAM_desc_t *descZ, CHAM_desc_t *descsum, double* theta,
                                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t Z = *descZ;
    CHAM_desc_t sum = *descsum;
    struct starpu_codelet *cl = &cl_ng_loglike;


    for (m = 0; m < Z.mt; m++) {
        tempmm = m == Z.mt - 1 ? Z.m - m * Z.mb : Z.mb;

        m0 = m * Z.mb;

        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descZ, ChamRealDouble, m, 0),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descsum, ChamRealDouble, 0, 0),
                           STARPU_VALUE, &theta[0], sizeof(double),
                           STARPU_VALUE, &theta[1], sizeof(double),
                           STARPU_VALUE, &theta[2], sizeof(double),
                           STARPU_VALUE, &theta[3], sizeof(double),
                           STARPU_VALUE, &theta[4], sizeof(double),
                           STARPU_VALUE, &theta[5], sizeof(double),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                           STARPU_NAME, "ng_loglike",
#endif
                           0);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    CHAMELEON_Sequence_Wait(sequence);
    return CHAMELEON_SUCCESS;
}