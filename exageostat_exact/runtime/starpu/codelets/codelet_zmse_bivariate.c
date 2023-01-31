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
 * StarPU codelet to Calculate Mean Square Error (MSE) between two vectors.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_dmse_bivariate_starpu(void *buffers[], void *cl_arg) {
    int m, m0, i;
    double* zpre;
    double* zmiss;
    double* serror1;
    double* serror2;
    double* serror;
    double local_serror1 = 0.0;
    double local_serror2 = 0.0;
    double local_serror = 0.0;

    serror1 = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    serror2 = (double* ) STARPU_MATRIX_GET_PTR(buffers[1]);
    serror = (double* ) STARPU_MATRIX_GET_PTR(buffers[2]);
    zpre = (double* ) STARPU_MATRIX_GET_PTR(buffers[3]);
    zmiss = (double* ) STARPU_MATRIX_GET_PTR(buffers[4]);

    starpu_codelet_unpack_args(cl_arg, &m, &m0);

    for (i = 0; i < m; i++) {
        if (i % 2 == 0) {
            local_serror1 += pow((zpre[i] - zmiss[i]), 2);
        } else
            local_serror2 += pow((zpre[i] - zmiss[i]), 2);
        local_serror += pow((zpre[i] - zmiss[i]), 2);
    }
    *serror1 += local_serror1;
    *serror2 += local_serror2;
    *serror += local_serror;
}

static struct starpu_codelet cl_dmse_bivariate =
        {
                .where        = STARPU_CPU,
                .cpu_funcs    = {CORE_dmse_bivariate_starpu},
                .nbuffers    = 5,
                .modes        = {STARPU_RW, STARPU_RW, STARPU_RW, STARPU_R, STARPU_R},
                .name        = "dmse_bivariate"
        };


/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_dmse_bivariate_Tile_Async - Calculate mean square eror (MSE) scalar value the prediction.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descZpre
 *           descZpre:  Observed measurements descZpre
 *
 * @param[in] descZmiss
 *           descZpre: Missing measurements descZpre
 *
 * @param[out] descserror
 *           descerror:  Mean Square Error (MSE)
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

int EXAGEOSTAT_MLE_dmse_bivariate_Tile_Async(CHAM_desc_t *descZpre, CHAM_desc_t *descZmiss, CHAM_desc_t *descserror1,
                                            CHAM_desc_t *descserror2, CHAM_desc_t *descserror,
                                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t Zpre = *descZpre;
    struct starpu_codelet *cl = &cl_dmse_bivariate;


    for (m = 0; m < Zpre.mt; m++) {
        tempmm = m == Zpre.mt - 1 ? Zpre.m - m * Zpre.mb : Zpre.mb;
        m0 = m * Zpre.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror1, ChamRealDouble, 0, 0),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror2, ChamRealDouble, 0, 0),
                           STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror, ChamRealDouble, 0, 0),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descZpre, ChamRealDouble, m, 0),
                           STARPU_R, EXAGEOSTAT_RTBLKADDR(descZmiss, ChamRealDouble, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                           STARPU_NAME, "dmse_bivariate",
#endif
                           0);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    CHAMELEON_Sequence_Wait(sequence);
    return CHAMELEON_SUCCESS;
}