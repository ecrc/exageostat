/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_ng_transform_lr.c
 *
 * StarPU codelet to Transform the measurements vector inside the non-Gaussian MLE function.
 *
 * @version 1.2.0
 *
 * @author Sagnik Mondal
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"
#include <hicma_struct.h>
#include <hicma_constants.h>
#include "hicma/hicma_ext/control/hicma_context.h"

static void CORE_ng_transform_lr_starpu(void *buffers[], void *cl_arg) {
    int m, m0;
    double* z;
    double* theta;
    double flag = 0;
    double* all_flags = &flag;
    *all_flags = 0;

    theta = (double* ) malloc(6 * sizeof(double));
    z = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    all_flags = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);   //to be removed

    starpu_codelet_unpack_args(cl_arg, &m, &m0, /*&flag,*/
                               &theta[0], &theta[1], &theta[2],
                               &theta[3], &theta[4], &theta[5]);

    double local_flag = 0.0;// = core_dmdet(A, m, n, m0, n0);

    //Transform the measurements vector inside the non-Gaussian MLE function.
    core_ng_transform(z, &local_flag, theta, m);
}

static struct starpu_codelet cl_ng_transform_lr =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_ng_transform_lr_starpu},
                .nbuffers    = 2,
                .modes        = {STARPU_RW, STARPU_W},
                .name        = "ng_transform_lr"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_ng_transform_lr_Tile_Async - Transform the measurements vector inside the non-Gaussian MLE function.
 *  Operates on vectors stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descZ
 *           descZ:  Observed measurements descZpre
 *
 * @param[in] theta
 *           theta: Model parameters.
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

int EXAGEOSTAT_ng_transform_lr_Tile_Async(HICMA_desc_t *descZ, HICMA_desc_t *descflag, const double* theta,
                                      HICMA_sequence_t *sequence, HICMA_request_t *request) {
    HICMA_context_t *hicmatxt;
    HICMA_option_t options;
    hicmatxt = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
	    return -2;
    HICMA_RUNTIME_options_init(&options, hicmatxt, sequence, request);
    int m, m0;
    int tempmm;
    HICMA_desc_t Z = *descZ;
    struct starpu_codelet *cl = &cl_ng_transform_lr;


    for (m = 0; m < Z.mt; m++) {
	    tempmm = m == Z.mt - 1 ? Z.m - m * Z.mb : Z.mb;
	    m0 = m * Z.mb;
	    starpu_insert_task(starpu_mpi_codelet(cl),
			    STARPU_VALUE, &tempmm, sizeof(int),
			    STARPU_VALUE, &m0, sizeof(int),
			    STARPU_RW, EXAGEOSTAT_HICMA_RTBLKADDR(descZ, HicmaRealDouble, m, 0),
			    STARPU_W, EXAGEOSTAT_HICMA_RTBLKADDR(descflag, HicmaRealDouble, 0, 0),
			    STARPU_VALUE, &theta[0], sizeof(double),
			    STARPU_VALUE, &theta[1], sizeof(double),
			    STARPU_VALUE, &theta[2], sizeof(double),
			    STARPU_VALUE, &theta[3], sizeof(double),
			    STARPU_VALUE, &theta[4], sizeof(double),
			    STARPU_VALUE, &theta[5], sizeof(double),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
			    STARPU_NAME, "ng_transform_lr",
#endif
			    0);
    }
    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, hicmatxt);
    //    HICMA_Sequence_Wait(sequence);
    return HICMA_SUCCESS;
}
