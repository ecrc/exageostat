/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dmat_reg.c
 *
 * StarPU codelet to Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2019-02-23
 *
 **/
#include "../include/starpu_exageostat.h"

static void cl_sdmat_reg_cpu_func(void *buffers[], void *cl_arg) {
    double* A;
    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
}

static struct starpu_codelet cl_sdmat_reg =
        {
                .where        = STARPU_CPU,
                .cpu_func    = cl_sdmat_reg_cpu_func,
                .nbuffers    = 1,
                .modes        = {STARPU_W},
                .name        = "dmat_reg"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_sdmat_reg_Tile_Async - Codelet to register the lower part of the matrix as souble and the upper one as float.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through buffersiptors.
 *  All dimensions are taken from the buffersiptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *		Upper or lower fill of the matrix.	    
 * 
 * @param[out] descA 
 *		descA:  Chameleon buffersiptor that handles the generated covariance matrix.
 *
 * @param[in] sequence
 *		Identifies the sequence of function calls that this call belongs to
 *		(for completion checks and exception handling purposes).
 *
 * @param[out] request
 *		Identifies this function call (for exception handling purposes).
 *
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
void EXAGEOSTAT_MLE_sdmat_reg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, RUNTIME_sequence_t *sequence,
                                        RUNTIME_request_t *request) {

    int n, m;

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_sdmat_reg;
    for (n = 0; n < A.nt; n++) {
        for (m = 0; m < A.mt; m++) {
            if (n > m)
                starpu_insert_task(starpu_mpi_codelet(cl),
                                   STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealFloat, m, n), 0);
            else
                starpu_insert_task(starpu_mpi_codelet(cl),
                                   STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, n), 0);
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}