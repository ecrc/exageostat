/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_zcorr.c
 *
 * StarPU codelet to calculate the r correlation vector (non-Gaussian kernel) 
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2021-03-04
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_dcorr_starpu(void *buffers[], void *cl_arg) {
    int i, m, m0;
    location *lmiss;
    location *lobs;
    double* theta;
    double* r;
    double expr = 0.0, x1, y1, x2, y2;
    double con = 0.0;
    double sigma_square = 1;//localtheta[0];// * localtheta[0];

    r = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &m0, &lmiss, &lobs, &theta);

    con = pow(2, (theta[1] - 1)) * tgamma(theta[1]);
    con = 1.0 / con;
    con = sigma_square * con;
    x1 = lmiss->x[0];
    y1 = lmiss->y[0];

    for (i = 0; i < m; i++) {
        x2 = lobs->x[m0 + i];
        y2 = lobs->y[m0 + i];
        expr = 4 * sqrt(2 * theta[1]) * (sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2))) / theta[0];

        if (expr == 0)
            r[i] = sigma_square /*+ 1e-4*/;
        else {
            r[i] = con * pow(expr, theta[1])
                   * gsl_sf_bessel_Knu(theta[1], expr); // Matern Function
        }
    }
}

static struct starpu_codelet cl_dcorr =
        {
                .where          = STARPU_CPU,
                .cpu_funcs      = {CORE_dcorr_starpu},
                .nbuffers       = 1,
                .modes          = {STARPU_W},
                .name           = "dcorr"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_MLE_dcorr_Async- codelet to compute dot product of A.A.
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
int EXAGEOSTAT_ng_corr_vec_gen_Tile_Async(
        CHAM_desc_t *descr, location *l, location *lobs,
        double* theta, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, m0;
    int tempmm;
    CHAM_desc_t r = *descr;
    struct starpu_codelet *cl = &cl_dcorr;

    for (m = 0; m < r.mt; m++) {
        tempmm = m == r.mt - 1 ? r.m - m * r.mb : r.mb;
        m0 = m * r.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_W, EXAGEOSTAT_RTBLKADDR(descr, ChamRealDouble, m, 0),
                           STARPU_VALUE, &l, sizeof(location *),
                           STARPU_VALUE, &lobs, sizeof(location *),
                           STARPU_VALUE, &theta, sizeof(double* ),
                           0);

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}