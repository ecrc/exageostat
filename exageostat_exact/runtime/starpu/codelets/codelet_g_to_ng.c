/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_g_to_ng.c
 *
 * StarPU codelet to Convert Gaussian measurements to non-Gaussian measurements.
 *
 * @version 1.2.0
 *
 * @author Sagnik Mondal
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_g_to_ng_starpu(void *buffers[], void *cl_arg) {
    int m, m0, i;
    double* z;
    double* theta;
    theta = (double* ) malloc(6 * sizeof(double));
    z = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &m0,
                               &theta[0], &theta[1], &theta[2],
                               &theta[3], &theta[4], &theta[5]);

    //core function to convert Z tile from Gaussian to non-Gaussian.
    core_g_to_ng(z, theta, m);
}

static struct starpu_codelet cl_g_to_ng =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_g_to_ng_starpu},
                .nbuffers    = 1,
                .modes        = {STARPU_RW},
                .name        = "g_to_ng"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_g_to_ng_Tile_Async - Convert Gaussian measurements to non-Gaussian measurements.
 *  Operates on Z vector stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descZ
 *           descZ:  Observed measurements descZ
 *
 * @param[in] theta
 *           theta: Model paramters
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

int EXAGEOSTAT_g_to_ng_Tile_Async(CHAM_desc_t *descZ, double* theta, RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t Z = *descZ;
    struct starpu_codelet *cl = &cl_g_to_ng;


    for (m = 0; m < Z.mt; m++) {
        tempmm = m == Z.mt - 1 ? Z.m - m * Z.mb : Z.mb;

        m0 = m * Z.mb;

        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descZ, ChamRealDouble, m, 0),
                           STARPU_VALUE, &theta[0], sizeof(double),
                           STARPU_VALUE, &theta[1], sizeof(double),
                           STARPU_VALUE, &theta[2], sizeof(double),
                           STARPU_VALUE, &theta[3], sizeof(double),
                           STARPU_VALUE, &theta[4], sizeof(double),
                           STARPU_VALUE, &theta[5], sizeof(double),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                           STARPU_NAME, "g_to_ng",
#endif
                           0);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    CHAMELEON_Sequence_Wait(sequence);

    return CHAMELEON_SUCCESS;
}
