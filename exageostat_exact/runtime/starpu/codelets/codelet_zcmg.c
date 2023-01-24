/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dcmg.c
 *
 * StarPU codelet to Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/starpu_exageostat.h"

static void cl_dcmg_cpu_func(void *buffers[], void *cl_arg) {
    int m, n, m0, n0;
    location *l1;
    location *l2;
    location *lm;
    double sigma1;
    double sigma2;
    double beta;
    double nu1;
    double nu2;
    double nu12;
    double a;
    double rho;
    double* theta;
    double* A;
    int distance_metric;
    int kernel;
    A = (double* ) STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &l1, &l2, &lm, &theta, &distance_metric, &kernel);//, &size);

    if (kernel == 0)
        core_dcmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 1)
        core_dcmg_nono_stat(A, m, n, m0, n0, l1, l2, lm, theta, distance_metric);
    else if (kernel == 2)
        core_dcmg_bivariate_flexible(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 3) //parsimonious
        core_dcmg_bivariate_parsimonious(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 4)
        core_dcmg_nuggets(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 5)
        core_dcmg_spacetime_matern(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 6)
        core_dcmg_matern_dsigma_square(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 7)
        core_dcmg_matern_dnu(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 8)
        core_dcmg_matern_dbeta(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 9)
        core_dcmg_matern_ddsigma_square(A, m, n);
    else if (kernel == 10)
        core_dcmg_matern_ddsigma_square_beta(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 11)
        core_dcmg_matern_ddsigma_square_nu(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 12)
        core_dcmg_matern_ddbeta_beta(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 13)
        core_dcmg_matern_ddbeta_nu(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 14)
        core_dcmg_matern_ddnu_nu(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 15)
        core_dcmg_spacetime_bivariate_parsimonious(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 16)
        core_ng_dcmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 17)
        core_ng_exp_dcmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 18)
        core_dcmg_spacetime_bivariate_parsimonious(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 19)
        core_dcmg_trivariate_parsimonious(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    else if (kernel == 20)
        core_dcmg_non_stat(A, m, n, m0, n0, l1, l2, theta, distance_metric);
}


#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_dcmg_cuda_func(void *buffers[], void *cl_arg)
{
    // Avoid generating matrix on GPU
    /*int m, n, m0, n0;
    location *l1;
    location *l2;
    double* theta;
    double* A;
    int distance_metric;
    theta   = (double* ) malloc(3* sizeof(double));
    A       = (double* )STARPU_MATRIX_GET_PTR(buffers[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &l1, &l2, &theta[0], &theta[1], &theta[2], &distance_metric);
    cudaStream_t stream = starpu_cuda_get_local_stream();
    cuda_dcmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
    cudaStreamSynchronize( stream );
*/
}
#endif


static struct starpu_codelet cl_dcmg =
        {
                .where        = STARPU_CPU /*| STARPU_CUDA*/,
                .cpu_func     = cl_dcmg_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)
                //    .cuda_func      = {cl_dcmg_cuda_func},
#endif
                .nbuffers     = 1,
                .modes        = {STARPU_W},
                .name         = "dcmg"
        };

//******************************************************************************
static void cl_scmg_cpu_func(void *buffers[], void *cl_arg) {
    int m, n, m0, n0;
    location *l1;
    location *l2;
    location *lm;
    double* theta;
    float *A;
    int distance_metric;
    int kernel;
    A = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &l1, &l2, &lm, &theta, &distance_metric, &kernel);

    if (kernel == 0)
        core_scmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
}


#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_scmg_cuda_func(void *buffers[], void *cl_arg)
{
int m, n, m0, n0;
location *l1;
location *l2;
double* theta;
double* A;
int distance_metric;
theta   = (double* ) malloc(3* sizeof(double));
A       = (double* )STARPU_MATRIX_GET_PTR(buffers[0]);
starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &l1, &l2, &theta[0], &theta[1], &theta[2], &distance_metric);
cudaStream_t stream = starpu_cuda_get_local_stream();
cuda_scmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
cudaStreamSynchronize( stream );
}
#endif


static struct starpu_codelet cl_scmg =
        {
                .where          = STARPU_CPU | STARPU_CUDA,
                .cpu_func       = cl_scmg_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)
                //    .cuda_func      = {cl_scmg_cuda_func},
#endif
                .nbuffers       = 1,
                .modes          = STARPU_W,
                .name           = "scmg"
        };

//******************************************************************************
static void cl_sdcmg_cpu_func(void *buffers[], void *cl_arg) {
    int m, n, m0, n0;
    location *l1;
    location *l2;
    location *lm;
    double* theta;
    float *A;
    int distance_metric;
    int kernel;
    A = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &l1, &l2, &lm, &theta, &distance_metric, &kernel);

    if (kernel == 0)
        core_sdcmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);
}

static struct starpu_codelet cl_sdcmg =
        {
                .where          = STARPU_CPU | STARPU_CUDA,
                .cpu_func       = cl_sdcmg_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)
                //    .cuda_func      = {cl_sdcmg_cuda_func},
#endif
                .nbuffers       = 1,
                .modes          = {STARPU_W},
                .name           = "sdcmg"
        };

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  EXAGEOSTAT_MLE_dcmg_Tile_Async - Codelet to generate covariance matrix in buffersiptor descA in  dense format between two sets of locations (l1, l2) (Matern Kernel).
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through buffersiptors.
 *  All dimensions are taken from the buffersiptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *        Upper or lower fill of the matrix.        
 * 
 * @param[out] descA 
 *        descA:  Chameleon buffersiptor that handles the generated covariance matrix.
 *
 * @param[in] sequence
 *        Identifies the sequence of function calls that this call belongs to
 *        (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *        Identifies this function call (for exception handling purposes).
 *
 * @param[in] l1
 *        Location struct of the first input.
 *
 * @param[in] l2
 *        Location struct of the second input.
 *
 * @param[in] theta
 *        Parameter vector that should be used to generate the output covariance matrix.
 *
 * @param[in] dm
 *        Distance metric "euclidean Distance ("ED" -->0) or "Great Circle Distance (GCD) -->1".
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
int EXAGEOSTAT_MLE_dcmg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                  location *l2, location *lm, double* theta,
                                  char* dm, char* kernel_fun,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, n, m0, n0;
    int distance_metric = strcmp(dm, "gc") == 0 ? 1 : 0;
    int kernel;
    if (strcmp(kernel_fun, "univariate_matern_stationary") == 0)
        kernel = 0;
    else if (strcmp(kernel_fun, "univariate_matern_non_stationary") == 0)
        kernel = 1;
    else if (strcmp(kernel_fun, "bivariate_matern_flexible") == 0)
        kernel = 2;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious") == 0
             || strcmp(kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        kernel = 3;
    else if (strcmp(kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        kernel = 4;
    else if (strcmp(kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        kernel = 5;
    else if (strcmp(kernel_fun, "univariate_matern_dsigma_square") == 0)
        kernel = 6;
    else if (strcmp(kernel_fun, "univariate_matern_dnu") == 0)
        kernel = 7;
    else if (strcmp(kernel_fun, "univariate_matern_dbeta") == 0)
        kernel = 8;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square") == 0)
        kernel = 9;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square_beta") == 0)
        kernel = 10;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square_nu") == 0)
        kernel = 11;
    else if (strcmp(kernel_fun, "univariate_matern_ddbeta_beta") == 0)
        kernel = 12;
    else if (strcmp(kernel_fun, "univariate_matern_ddbeta_nu") == 0)
        kernel = 13;
    else if (strcmp(kernel_fun, "univariate_matern_ddnu_nu") == 0)
        kernel = 14;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        kernel = 15;
    else if (strcmp(kernel_fun, "univariate_matern_non_gaussian") == 0)
        kernel = 16;
    else if (strcmp(kernel_fun, "univariate_exp_non_gaussian") == 0)
        kernel = 17;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        kernel = 18;
    else if (strcmp(kernel_fun, "trivariate_matern_parsimonious") == 0 ||
             strcmp(kernel_fun, "trivariate_matern_parsimonious_profile") == 0)
        kernel = 19;
    else if (strcmp(kernel_fun, "univariate_matern_non_stat") == 0)
        kernel = 20;
    else {
        fprintf(stderr, "Choosen kernel is not exist: %s!\n", kernel_fun);
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    int tempmm, tempnn;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_dcmg;
    int size = A.n;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (uplo == ChamUpperLower)
            m = 0;
        else
            m = A.m == A.n ? n : 0;
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, n),
                               STARPU_VALUE, &l1, sizeof(location * ),
                               STARPU_VALUE, &l2, sizeof(location * ),
                               STARPU_VALUE, &lm, sizeof(location * ),
                               STARPU_VALUE, &theta, sizeof(double* ),
                               STARPU_VALUE, &distance_metric, sizeof(int),
                               STARPU_VALUE, &kernel, sizeof(int),
                               0);
        }

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}

/*******************************************************************************
 *
 * @ingroup CHAMELEON_Complex32_t_Tile
 *
 *  EXAGEOSTAT_MLE_scmg_Tile_Async - Calculate covariance matrix descA.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through buffersiptors.
 *  All dimensions are taken from the buffersiptors.
 *
 *******************************************************************************
 *
 * @param[out] descA
 *           descA:  Chameleon buffersiptor that handles the generated covariance matrix
 *
 * @param[in] l1
 *           location struct of the first input
 *
 * @param[in] l2
 *          location struct of the second input
 *
 * @param[in] theta
 *           parameter vector that should be used to generate the output covariance matrix
 *
 * @param[in] dm
 *           distance metric "euclidean Distance (ED"" or "Great Circle Distance (GCD)"
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


int EXAGEOSTAT_MLE_scmg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                  location *l2, location *lm, double* theta,
                                  char* dm, char* kernel_fun,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, n, m0, n0;
    int distance_metric = strcmp(dm, "gc") == 0 ? 1 : 0;
    int kernel;
    if (strcmp(kernel_fun, "univariate_matern_stationary") == 0)
        kernel = 0;
    else if (strcmp(kernel_fun, "univariate_matern_non_stationary") == 0)
        kernel = 1;
    else if (strcmp(kernel_fun, "bivariate_matern_flexible") == 0)
        kernel = 2;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious") == 0
             || strcmp(kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        kernel = 3;
    else if (strcmp(kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        kernel = 4;
    else if (strcmp(kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        kernel = 5;
    else if (strcmp(kernel_fun, "univariate_matern_dsigma_square") == 0)
        kernel = 6;
    else if (strcmp(kernel_fun, "univariate_matern_dnu") == 0)
        kernel = 7;
    else if (strcmp(kernel_fun, "univariate_matern_dbeta") == 0)
        kernel = 8;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square") == 0)
        kernel = 9;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square_beta") == 0)
        kernel = 10;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square_nu") == 0)
        kernel = 11;
    else if (strcmp(kernel_fun, "univariate_matern_ddbeta_beta") == 0)
        kernel = 12;
    else if (strcmp(kernel_fun, "univariate_matern_ddbeta_nu") == 0)
        kernel = 13;
    else if (strcmp(kernel_fun, "univariate_matern_ddnu_nu") == 0)
        kernel = 14;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        kernel = 15;
    else if (strcmp(kernel_fun, "univariate_matern_non_gaussian") == 0)
        kernel = 16;
    else if (strcmp(kernel_fun, "univariate_exp_non_gaussian") == 0)
        kernel = 17;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        kernel = 18;
    else if (strcmp(kernel_fun, "trivariate_matern_parsimonious") == 0 ||
             strcmp(kernel_fun, "trivariate_matern_parsimonious_profile") == 0)
        kernel = 19;
    else if (strcmp(kernel_fun, "univariate_matern_non_stat") == 0)
        kernel = 20;

    else {
        fprintf(stderr, "Choosen kernel is not exist(5)!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    int tempmm, tempnn;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *cl = &cl_scmg;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (uplo == ChamUpperLower)
            m = 0;
        else
            m = A.m == A.n ? n : 0;
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealFloat, m, n),
                               STARPU_VALUE, &l1, sizeof(location * ),
                               STARPU_VALUE, &l2, sizeof(location * ),
                               STARPU_VALUE, &lm, sizeof(location * ),
                               STARPU_VALUE, &theta, sizeof(double* ),
                               STARPU_VALUE, &distance_metric, sizeof(int),
                               STARPU_VALUE, &kernel, sizeof(int),
                    //STARPU_VALUE, &num_locs, sizeof(int),
                    //STARPU_VALUE, &num_params, sizeof(int),
                               0);
        }

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}

/*******************************************************************************/
int EXAGEOSTAT_MLE_sdcmg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                   location *l2, location *lm, double* theta,
                                   char* dm, char* kernel_fun,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    int m, n, m0, n0;
    int distance_metric = strcmp(dm, "gc") == 0 ? 1 : 0;
    int kernel;

    if (strcmp(kernel_fun, "univariate_matern_stationary") == 0)
        kernel = 0;
    else if (strcmp(kernel_fun, "univariate_matern_non_stationary") == 0)
        kernel = 1;
    else if (strcmp(kernel_fun, "bivariate_matern_flexible") == 0)
        kernel = 2;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious") == 0
             || strcmp(kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        kernel = 3;
    else if (strcmp(kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        kernel = 4;
    else if (strcmp(kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        kernel = 5;
    else if (strcmp(kernel_fun, "univariate_matern_dsigma_square") == 0)
        kernel = 6;
    else if (strcmp(kernel_fun, "univariate_matern_dnu") == 0)
        kernel = 7;
    else if (strcmp(kernel_fun, "univariate_matern_dbeta") == 0)
        kernel = 8;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square") == 0)
        kernel = 9;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square_beta") == 0)
        kernel = 10;
    else if (strcmp(kernel_fun, "univariate_matern_ddsigma_square_nu") == 0)
        kernel = 11;
    else if (strcmp(kernel_fun, "univariate_matern_ddbeta_beta") == 0)
        kernel = 12;
    else if (strcmp(kernel_fun, "univariate_matern_ddbeta_nu") == 0)
        kernel = 13;
    else if (strcmp(kernel_fun, "univariate_matern_ddnu_nu") == 0)
        kernel = 14;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        kernel = 15;
    else if (strcmp(kernel_fun, "univariate_matern_non_gaussian") == 0)
        kernel = 16;
    else if (strcmp(kernel_fun, "univariate_exp_non_gaussian") == 0)
        kernel = 17;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        kernel = 18;
    else if (strcmp(kernel_fun, "trivariate_matern_parsimonious") == 0 ||
             strcmp(kernel_fun, "trivariate_matern_parsimonious_profile") == 0)
        kernel = 19;
    else if (strcmp(kernel_fun, "univariate_matern_non_stat") == 0)
        kernel = 20;
    else {
        fprintf(stderr, "Choosen kernel is not exist(3)!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }


    int tempmm, tempnn;
    CHAM_desc_t A = *descA;
    struct starpu_codelet *dcl = &cl_dcmg;
    struct starpu_codelet *sdcl = &cl_sdcmg;

    int k = 0;
    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (uplo == ChamUpperLower)
            m = 0;
        else
            m = A.m == A.n ? n : 0;
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(dcl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, n),
                               STARPU_VALUE, &l1, sizeof(location * ),
                               STARPU_VALUE, &l2, sizeof(location * ),
                               STARPU_VALUE, &lm, sizeof(location * ),
                               STARPU_VALUE, &theta, sizeof(double* ),
                               STARPU_VALUE, &distance_metric, sizeof(int),
                               STARPU_VALUE, &kernel, sizeof(int),
                               0);
        }

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}
