/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *  @file MLE_lr.c
 *
 *  
 *  ExaGeoStat is a software package provided by KAUST,
 *  King Abdullah University of Science and Technology - ECRC
 *
 * @version 1.2.0
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/MLE_lr.h"
#include "hicma/misc/auxcompute_z.h"
#include "hicma_d.h"

//***************************************************************************************
//HiCMA global variables.
int store_only_diagonal_tiles = 1;
STARSH_blrf *mpiF;
int use_scratch = 1;
int global_check = 0;  //used to create dense matrix for accuracy check
int calc_rank_stat = 0; //Calculate the rank. Set to 0 in normal execution.
double *Ark_initial;
double *Ark_old;

int
EXAGEOSTAT_TLR_MLE_dzvg_Tile(MLE_data *data, double *Nrand, double *initial_theta, int n, int dts, int log, int p_grid,
                             int q_grid) {
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- HICMA-sync
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] dts: tile size (MB) is used only in the case of HiCMA not CHAM.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * @param[in] p_grid: p_grid in the case of distributed system.
     * @param[in] q_grid: q_grid in the case of distributed system.
     * */
    //Initialization
    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAMELEON_descZ = NULL;
    CHAM_desc_t *CHAMELEON_descproduct = NULL;    
    double *univariate_theta;
    double *univariate2_theta;
    double *univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    //Create two dense descriptors to generate the measurment vectors
    CHAMELEON_Sequence_Create(&msequence);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC, NULL, ChamRealDouble, dts, dts, dts * dts, n, n, 0, 0, n, n, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZ, NULL, ChamRealDouble, dts, dts, dts * dts, n, 1, 0, 0, n, 1,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descproduct, &data->dotp, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1,
                                    1, p_grid, q_grid);

    //Generate the co-variance matrix C
    VERBOSE("LR: Initializing Covariance Matrix (Synthetic Dataset Dense Generation Phase).....");

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {

        univariate_theta = (double *) malloc(3 * sizeof(double));
        univariate2_theta = (double *) malloc(3 * sizeof(double));
        univariate3_theta = (double *) malloc(3 * sizeof(double));
        univariate_theta[0] = initial_theta[0];
        univariate_theta[1] = initial_theta[2];
        univariate_theta[2] = initial_theta[3];

        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, chameleon_desc_submatrix(CHAM_descC, 0, 0, CHAM_descC->m / 2,
                                                                                CHAM_descC->n / 2), &data->l1,
                                       &data->l1,
                                       &data->lm, univariate_theta, data->dm,
                                       "univariate_matern_stationary", msequence, mrequest);


        nu12 = 0.5 * (initial_theta[3] + initial_theta[4]);
        rho = initial_theta[5] * sqrt((tgamma(initial_theta[3] + 1) * tgamma(initial_theta[4] + 1)) /
                                      (tgamma(initial_theta[3]) * tgamma(initial_theta[4]))) *
              tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(initial_theta[0] * initial_theta[1]);

        univariate2_theta[0] = sigma_square12;
        univariate2_theta[1] = initial_theta[2];
        univariate2_theta[2] = nu12;
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, chameleon_desc_submatrix(CHAM_descC, CHAM_descC->m / 2,
                                                                                0, CHAM_descC->m / 2,
                                                                                CHAM_descC->n / 2),
                                       &data->l1, &data->l1, &data->lm,
                                       univariate2_theta, data->dm, "univariate_matern_stationary", msequence,
                                       mrequest);

        univariate3_theta[0] = initial_theta[1];
        univariate3_theta[1] = initial_theta[2];
        univariate3_theta[2] = initial_theta[4];
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, chameleon_desc_submatrix(CHAM_descC, CHAM_descC->m / 2,
                                                                                CHAM_descC->n / 2, CHAM_descC->m / 2,
                                                                                CHAM_descC->n / 2), &data->l1,
                                       &data->l1,
                                       &data->lm, univariate3_theta, data->dm, "univariate_matern_stationary",
                                       msequence,
                                       mrequest);
    } else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0) {
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1,
                                       &data->l1, &data->lm, initial_theta,
                                       data->dm, "univariate_matern_stationary", msequence, mrequest);
    } else {
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC, &data->l1,
                                       &data->l1, &data->lm, initial_theta,
                                       data->dm, data->kernel_fun, msequence, mrequest);

    }

    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("LR: Generate Normal Random Distribution Vector Z (Synthetic Dataset Dense Generation Phase) .....");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAMELEON_descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("LR: Cholesky factorization of Sigma (Synthetic Dataset Dense Generation Phase) .....");
    int success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_descC);
    SUCCESS(success, "LR: Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("LR: Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Dense Generation Phase) .....");
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAMELEON_descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if (log == 1) {
        double *z;
#if defined(CHAMELEON_USE_MPI)
        z = (double* ) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAMELEON_descZ, z, n);
#else
        z = CHAMELEON_descZ->mat;
#endif
        write_vectors(z, data, n);
#if defined(CHAMELEON_USE_MPI)
        free(z);
#endif
    }

    VERBOSE("LR: Done Z Vector Dense Generation Phase. (Chameleon Synchronous)\n");
    VERBOSE("************************************************************\n");

    data->descZ = CHAMELEON_descZ;
    data->descproduct = CHAMELEON_descproduct;
    data->sequence = msequence;
    data->request = mrequest;
    //Destory dense descriptors to save memory.
    CHAMELEON_Desc_Destroy(&CHAM_descC);
    return 0;
}


void EXAGEOSTAT_TLR_MLE_zcpy(MLE_data *data, double *streamdata)
//! Copy measurements vector from Lapack
/*! format to Chameleon format.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] streamdata: measurments vector in lapack format.
 * */
{
    RUNTIME_sequence_t *hsequence = (RUNTIME_sequence_t *) data->hsequence;
    RUNTIME_request_t *hrequest = (RUNTIME_request_t *) data->hrequest;
    VERBOSE("LR: Copy Z from vector to decriptor.\n");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(data->hicma_descZ, streamdata, hsequence, hrequest);
    CHAMELEON_Sequence_Wait(hsequence);
    VERBOSE("LR: Done Z copying step.\n");
    VERBOSE("************************************************************\n");
}

int EXAGEOSTAT_TLR_MLE_dzvg_Tile_Async(MLE_data *data, double *Nrand, double *initial_theta, int n, int dts, int log,
                                       int p_grid,
                                       int q_grid) {
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- HICMA-async
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] dts: tile size (MB) is used only in the case of HiCMA not CHAM.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * @param[in] p_grid: p_grid in the case of distributed system.
     * @param[in] q_grid: q_grid in the case of distributed system.
     * */
    //Initialization
    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAMELEON_descZ = NULL;

    //Create two dense descriptors to generate the measurment vectors
    CHAMELEON_Sequence_Create(&msequence);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC, NULL, ChamRealDouble, dts, dts, dts * dts, n, n, 0, 0, n, n, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZ, NULL, ChamRealDouble, dts, dts, dts * dts, n, 1, 0, 0, n, 1,
                                    p_grid, q_grid);

    //Generate the co-variance matrix C
    VERBOSE("LR: Initializing Covariance Matrix (Synthetic Dataset Dense Generation Phase).....");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm,
                                   data->kernel_fun, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("LR: Generate Normal Random Distribution Vector Z (Synthetic Dataset Dense Generation Phase) .....");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAMELEON_descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("LR: Cholesky factorization of Sigma (Synthetic Dataset Dense Generation Phase) .....");
    int success = CHAMELEON_dpotrf_Tile_Async(ChamLower, CHAM_descC, msequence, mrequest);
    SUCCESS(success, "LR: Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("LR: Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Dense Generation Phase) .....");
    CHAMELEON_dtrmm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAMELEON_descZ, msequence,
                               mrequest);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if (log == 1) {
        double *z;
#if defined(CHAMELEON_USE_MPI)
        z = (double* ) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAMELEON_descZ, z, n);
#else
        z = CHAMELEON_descZ->mat;
#endif
        write_vectors(z, data, n);
#if defined(CHAMELEON_USE_MPI)
        free(z);
#endif
    }

    VERBOSE("LR: Done Z Vector Dense Generation Phase. (Chameleon Synchronous)\n");
    VERBOSE("************************************************************\n");


    data->descZ = CHAMELEON_descZ;
    data->sequence = msequence;
    data->request = mrequest;
    //Destory dense descriptors to save memory.
    CHAMELEON_Desc_Destroy(&CHAM_descC);
    return 0;
}


//compute MLE function
double EXAGEOSTAT_TLR_dmle_Tile(unsigned n, const double *theta, double *grad, void *app_data)
//! Maximum Likelihood Evaluation (MLE)
/*!  -- HICMA-sync
 * Returns the loglikelihhod value for the given theta.
 * @param[in] n: unsigned variable used by NLOPT library.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] grad: double variable used by NLOPT library.
 * @param[in] app_data: MLE_data struct with different MLE inputs.
 * */
{
    //Initialization
    double loglik = 0.0;
    int compress_diag = 0;
    double logdet = 0.0;
    double time_facto = 0.0,
            time_solve = 0.0,
            logdet_calculate = 0.0,
            matrix_gen_time = 0.0,
            test_time = 0.0;
    double flops = 0.0;
    int success, maxrank, acc;
    int N, NRHS;
    int lts;//, dts;
    int hicma_data_type;
    int i = 0;
    int num_params = 0;
    MLE_data *data = (MLE_data *) app_data;
    HICMA_desc_t *hicma_descC = (HICMA_desc_t *) data->hicma_descC;
    HICMA_desc_t *hicma_descCD = (HICMA_desc_t *) data->hicma_descCD;
    HICMA_desc_t *hicma_descCUV = (HICMA_desc_t *) data->hicma_descCUV;
    HICMA_desc_t *hicma_descCrk = (HICMA_desc_t *) data->hicma_descCrk;
    HICMA_desc_t *hicma_descZ = (HICMA_desc_t *) data->hicma_descZ;
    CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) data->descZ;
    HICMA_desc_t *hicma_descZcpy = (HICMA_desc_t *) data->hicma_descZcpy;
    HICMA_desc_t *hicma_descdet = (HICMA_desc_t *) data->hicma_descdet;
    CHAM_desc_t *cham_descproduct = (CHAM_desc_t *) data->descproduct;
    HICMA_sequence_t *hsequence = (HICMA_sequence_t *) data->hsequence;
    HICMA_request_t *hrequest = (HICMA_request_t *) data->hrequest;
    if (strcmp(data->kernel_fun, "univariate_matern_stationary") == 0 ||
        strcmp(data->kernel_fun, "univariate_pow_exp_stationary") == 0)
        num_params = 3;
    else if (strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        num_params = 4;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0)
        num_params = 9;
    else if (strcmp(data->kernel_fun, "bivariate_matern_flexible") == 0)
        num_params = 11;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_gaussian") == 0 ||
             strcmp(data->kernel_fun, "univariate_exp_non_gaussian") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0)
        num_params = 10;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stat") == 0)
        num_params = 8;
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    //verbose               = data->verbose;
    N = hicma_descCD->m;
    NRHS = hicma_descZ->n;
    lts = hicma_descZ->mb;
    maxrank = ((MLE_data *) data)->hicma_maxrank;
    acc = ((MLE_data *) data)->hicma_acc;
    hicma_data_type = ((MLE_data *) data)->hicma_data_type;
    data->det = 0;
    data->dotp = 0;

    //Save a copy of descZ into descZcpy for restoring each iteration
    if (data->iter_count == 0) {

        if (strcmp(data->locsFPath, "") == 0) {
            double *z = (double *) malloc(N * sizeof(double));
            CHAMELEON_Tile_to_Lapack(CHAM_descZ, z, N);
            HICMA_Lapack_to_Tile(z, N, hicma_descZ);
            free(z);
            VERBOSE("transforming\n");
            //CHAMELEON_Desc_Destroy(&CHAM_descZ);
        }
        HICMA_dlacpy_Tile(ChamUpperLower, hicma_descZ, hicma_descZcpy);
    }
    if (strcmp(data->recovery_file, "") != 0 &&
        recover(data->recovery_file, data->iter_count, theta, &loglik, num_params));

    else {
        //Save a copy of descZ into descZcpy for restoring each iteration

        //Matrix generation part.       
        VERBOSE("LR:Generate New Covariance Matrix...");
        START_TIMING(matrix_gen_time);

        HICMA_problem_t hicma_problem;
        hicma_problem.theta = (double *) theta;
        hicma_problem.noise = 1e-4;
        hicma_problem.ndim = 2;

        // I need to change the location struct to vector array (TO DO)
        double *xycoord = (double *) malloc(2 * N * sizeof(double));

        if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0) {
            int j = 0;
            for (i = 0; i < N; i++) {
                xycoord[i] = data->l1.x[j];
                xycoord[N + i] = data->l1.y[j];
                if (i % 2 != 0)
                    j++;
            }
            hicma_problem.kernel_type =
                    strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_PARSIMONIOUS_GCD : STARSH_SPATIAL_PARSIMONIOUS_SIMD;
        } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0) {
            int j = 0;
            for (i = 0; i < N; i++) {
                xycoord[i] = data->l1.x[j];
                xycoord[N + i] = data->l1.y[j];
                if (i % 3 != 0)
                    j++;
            }
            hicma_problem.kernel_type =
                    strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_PARSIMONIOUS_GCD : STARSH_SPATIAL_PARSIMONIOUS_SIMD;
        } else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
                   strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
            int j = 0;
            for (i = 0; i < N; i++) {
                if (i == N / 2)
                    j = 0;
                xycoord[i] = data->l1.x[j];
                xycoord[N + i] = data->l1.y[j];
                j++;
            }
            hicma_problem.kernel_type =
                    strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_PARSIMONIOUS2_GCD : STARSH_SPATIAL_PARSIMONIOUS2_SIMD;
        } else {
            for (i = 0; i < N; i++) {
                xycoord[i] = data->l1.x[i];
                xycoord[N + i] = data->l1.y[i];
            }

            hicma_problem.kernel_type =
                    strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_MATERN2_GCD : STARSH_SPATIAL_MATERN2_SIMD;
        }
        hicma_problem.point = xycoord;
        HICMA_zgenerate_problem(hicma_data_type, 'S', 0, N, lts, hicma_descCUV->mt, hicma_descCUV->nt, &hicma_problem);
        mpiF = hicma_problem.starsh_format;

        HICMA_zgytlr_Tile(HicmaLower, hicma_descCUV, hicma_descCD, hicma_descCrk, 0, maxrank, pow(10, -1.0 * acc),
                          compress_diag, hicma_descC);
        STOP_TIMING(matrix_gen_time);
        VERBOSE(" Done.\n");
        fflush(stdout);
        //******************************
        VERBOSE("LR: re-Copy z...");
        START_TIMING(test_time);
        //re-store old Z
        HICMA_dlacpy_Tile(HicmaUpperLower, hicma_descZcpy, hicma_descZ);
        STOP_TIMING(test_time);
        VERBOSE(" Done.\n");
        //*************

        //Calculate Cholesky Factorization (C=LL-1)
        VERBOSE("LR: Cholesky factorization of Sigma...");
        START_TIMING(time_facto);
        success = HICMA_dpotrf_Tile(HicmaLower, hicma_descCUV, hicma_descCD, hicma_descCrk, 0, maxrank,
                                    pow(10, -1.0 * acc));
        SUCCESS(success, "LR: Factorization cannot be performed..\n The matrix is not positive definite\n\n");
        STOP_TIMING(time_facto);
        flops = flops + FLOPS_DPOTRF(N);
        VERBOSE(" Done.");

        //*********
        //Calculate log(|C|) --> log(square(|L|))
        VERBOSE("LR:Calculating the log determinant ...");
        START_TIMING(logdet_calculate);
        data->det = 0;
        EXAGEOSTAT_TLR_MLE_dmdet_Tile_Async(hicma_descCD, hsequence, &hrequest[0], hicma_descdet);
        HICMA_Sequence_Wait(hsequence);
        logdet = 2 * data->det;
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.");

        //Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("LR:Solving the linear system ...");
        START_TIMING(time_solve);
        //Compute triangular solve LC*X = Z

        HICMA_dtrsmd_Tile(HicmaLeft, HicmaLower, HicmaNoTrans, HicmaNonUnit, 1, hicma_descCUV, hicma_descCD,
                          hicma_descCrk,hicma_descZ, maxrank);
        STOP_TIMING(time_solve);
        flops = flops + FLOPS_DTRSM(ChamLeft, N, NRHS);
        VERBOSE(" Done.");
        VERBOSE("Copy to chameleon");

        double *z = (double *) malloc(N * N * sizeof(double));
        HICMA_Tile_to_Lapack(hicma_descZ, z, N);
        CHAMELEON_Lapack_to_Tile(z, N, CHAM_descZ);
        free(z);

        VERBOSE("LR:Calculating dot product...");
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, cham_descproduct);
        loglik = -0.5 * ((MLE_data *) data)->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);

        VERBOSE(" Done.");


        //Calculate the ranks
        if (calc_rank_stat == 1) {
            int i = 0;
            int MBrk = 1;
            int NBrk = 1;
            int Mrk = hicma_descCUV->mt;
            int Nrk = hicma_descCUV->mt;
            if (data->iter_count == 0) {
                Ark_initial = (double *) calloc(Mrk * Nrk, sizeof(double));
                Ark_old = (double *) calloc(Mrk * Nrk, sizeof(double));
            }
            HICMA_Tile_to_Lapack(hicma_descCrk, Ark_initial, Mrk);
            for (i = 0; i < Mrk * Nrk; i++) {
                if (Ark_initial[i] > Ark_old[i]) {
                    Ark_old[i] = Ark_initial[i];
                }
            }

            if (HICMA_My_Mpi_Rank() == 0) {
                fwrite_array(hicma_descCrk->m, hicma_descCrk->n, hicma_descCrk->m, Ark_old, "ranks.csv");
                print_array(hicma_descCrk->m, hicma_descCrk->n, hicma_descCrk->m, Ark_old, stdout);
                HICMA_stat_t hicma_statrk_initial;
                zget_stat(ChamLower, Ark_old, Mrk, Nrk, Mrk, &hicma_statrk_initial);
                zprint_stat(hicma_statrk_initial);
                fflush(stderr);
                fflush(stdout);
            }
        }

        //multiplicative scale
        data->variance = theta[0];
    }
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if(HICMA_My_Mpi_Rank() == 0)
    {
#endif

    printf(" %3d- Model Parameters (", ((MLE_data *) data)->iter_count + 1);

    if (data->log == 1)
        fprintf(((MLE_data *) data)->pFileLog, " %3d- Model Parameters (", ((MLE_data *) data)->iter_count + 1);

    if (strcmp(((MLE_data *) data)->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
        strcmp(((MLE_data *) data)->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        printf("%.8f, %.8f,", ((MLE_data *) data)->variance1, ((MLE_data *) data)->variance2);
        i = 2;

    } else
        i = 0;
    for (; i < num_params; i++) {
        printf("%.8f", theta[i]);
        if (i < num_params - 1)
            printf(",");

        if (((MLE_data *) data)->log == 1)
            fprintf(((MLE_data *) data)->pFileLog, "%.8f, ", theta[i]);
    }

    printf(")----> LogLi: %.18f\n", loglik);
    if (((MLE_data *) data)->log == 1)
        fprintf(((MLE_data *) data)->pFileLog, ")----> LogLi: %.18f\n", loglik);

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    ((MLE_data *) data)->iter_count++;
    // for experiments
    ((MLE_data *) data)->avg_exec_time_per_iter += matrix_gen_time + time_facto + logdet_calculate + time_solve;
    ((MLE_data *) data)->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    ((MLE_data *) data)->final_loglik = loglik;

    results.final_loglik = loglik;
    return loglik;
}

double EXAGEOSTAT_TLR_dmle_Tile_Async(unsigned n, const double *theta, double *grad, void *data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- HICMA-async
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
}

void
EXAGEOSTAT_TLR_dmle_Predict_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs, int lts, int p_grid, int q_grid,
                                     int mse_flag)
//! Allocate prediction operation descriptors.
/*!
 * Returns MLE_data data with initial values and new descriptors locations.
 * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] dts: tile size (MB).
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
 * */
{
    int MBC, NBC, MC, NC;
    int MBD, NBD, MD, ND;
    int MBUV, NBUV, MUV, NUV;
    int MBrk, NBrk, Mrk, Nrk;
    CHAM_desc_t *CHAMELEON_descZmiss = NULL;
    CHAM_desc_t *hicma_descC = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *hicma_descC22D = NULL;
    CHAM_desc_t *hicma_descC22UV = NULL;
    CHAM_desc_t *hicma_descC22rk = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAMELEON_descZactual = NULL;
    CHAM_desc_t *CHAMELEON_descZobs = NULL;

    MLE_data *data = (MLE_data *) CHAM_data;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return;
    }

    //Descriptors Creation
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZobs, NULL, ChamRealDouble, lts, lts, lts * lts, nZobs, 1, 0, 0,
                                    nZobs, 1, p_grid, q_grid);

    if (mse_flag == 1) {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZactual, NULL, ChamRealDouble, lts, lts, lts * lts, nZmiss, 1, 0,
                                        0, nZmiss, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealDouble, lts, lts, lts * lts, 1, 1, 0, 0,
                                        1, 1, p_grid, q_grid);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZmiss, NULL, ChamRealDouble, lts, lts, lts * lts, nZmiss, 1, 0, 0,
                                    nZmiss, 1, p_grid, q_grid);


    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC12, NULL, ChamRealDouble, lts, lts, lts * lts, nZmiss, nZobs, 0, 0,
                                    nZmiss, nZobs, p_grid, q_grid);
    //**********************************************************

    //Sameh
    //CDense Descriptor
    if (data->check == 1) {
        MBC = lts;
        NBC = lts;
        MC = nZobs;
        NC = nZobs;
    } else {
        MBC = 1;
        NBC = 1;
        MC = lts;
        NC = lts;
    }

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC, NULL, ChamRealDouble, MBC, NBC, MBC * NBC, MC, NC, 0, 0, MC, NC,
                                    p_grid, q_grid);
    MBD = lts;
    NBD = lts;
    MD = nZobs;
    ND = MBD;

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22D, NULL, ChamRealDouble, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND,
                                    p_grid, q_grid);
    //CAD Descriptor
    MBUV = lts;
    NBUV = 2 * data->hicma_maxrank;
    MUV = -1;
    int N_over_lts_times_lts = nZobs / lts * lts;
    if (N_over_lts_times_lts < nZobs) {
        MUV = N_over_lts_times_lts + lts;
    } else if (N_over_lts_times_lts == nZobs) {
        MUV = N_over_lts_times_lts;
    } else {
        printf("%s %d: This case should not happen\n", __FILE__, __LINE__);
        exit(-1);
    }
    double expr = (double) MUV / (double) lts;
    NUV = 2 * expr * data->hicma_maxrank;
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22UV, NULL, ChamRealDouble, MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0,
                                    MUV, NUV, p_grid, q_grid);

    //CUV Descriptor
    MBrk = 1;
    NBrk = 1;
    Mrk = hicma_descC22UV->mt;
    Nrk = hicma_descC22UV->mt;
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22rk, NULL, ChamRealDouble, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,
                                    Mrk, Nrk, p_grid, q_grid);

    //Initiate data descriptors
    data->descZmiss = CHAMELEON_descZmiss;
    data->hicma_descC = hicma_descC;
    data->hicma_descC22D = hicma_descC22D;
    data->hicma_descC22UV = hicma_descC22UV;
    data->hicma_descC22rk = hicma_descC22rk;
    data->descC12 = CHAM_descC12;
    data->descmse = CHAM_descmse;
    data->descZactual = CHAMELEON_descZactual;
    data->descZobs = CHAMELEON_descZobs;

}

//Predict missing values base on a set of given values and covariance matrix
double EXAGEOSTAT_TLR_dmle_Predict_Tile(MLE_data *CHAM_data, double *theta, int nZmiss, int nZobs, double *Zobs,
                                        double *Zactual,
                                        double *Zmiss, int n, int lts)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- HICMA-sync
 * Returns the prediction Mean Square Error (MSE) as double
 * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] Zobs: observed values vector (known observations).
 * @param[in] Zmiss missing values vector (unknown observations).
 * @param[in] Zactual: actual missing values vector (in the case of testing MSE).
 * @param[in] n: number of spatial locations.
 * */
{
    //initialization

    double time_solve = 0.0;
    double mat_gen_time = 0.0;
    double time_gemm = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;
    int maxrank = 0;
    int acc = 0;
    int hicma_data_type = 0;
    int compress_diag = 0;

    CHAM_desc_t *CHAMELEON_descZmiss = NULL;
    CHAM_desc_t *hicma_descC = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *hicma_descC22D = NULL;
    CHAM_desc_t *hicma_descC22UV = NULL;
    CHAM_desc_t *hicma_descC22rk = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAMELEON_descZactual = NULL;
    CHAM_desc_t *CHAMELEON_descZobs = NULL;

    MLE_data *data = (MLE_data *) CHAM_data;
    RUNTIME_sequence_t *hsequence = (RUNTIME_sequence_t *) data->hsequence;
    RUNTIME_request_t *hrequest = (RUNTIME_request_t *) data->hrequest;
    data->mserror = 0;
    maxrank = data->hicma_maxrank;
    acc = data->hicma_acc;
    hicma_data_type = data->hicma_data_type;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return -1;
    }

    //Initiate data descriptors
    CHAMELEON_descZmiss = data->descZmiss;
    hicma_descC = data->hicma_descC;
    hicma_descC22D = data->hicma_descC22D;
    hicma_descC22UV = data->hicma_descC22UV;
    hicma_descC22rk = data->hicma_descC22rk;
    CHAM_descC12 = data->descC12;
    CHAM_descmse = data->descmse;
    CHAMELEON_descZactual = data->descZactual;
    CHAMELEON_descZobs = data->descZobs;
    //****

    //Copy data to vectors
    VERBOSE("Copy measurments vector to descZobs descriptor...");

    CHAMELEON_Lapack_to_Tile(Zobs, nZobs, CHAMELEON_descZobs);
    VERBOSE(" Done.\n");

    if (Zactual != NULL) {
        VERBOSE("Copy actual measurments vector to descZactual descriptor...");
        CHAMELEON_Lapack_to_Tile(Zactual, nZmiss, CHAMELEON_descZactual);
        VERBOSE(" Done.\n");
    }
    //*********************************************

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
    theta[0] = data->variance;
    printf("estimated parameters: %f - %f - %f\n", theta[0], theta[1], theta[2]);
    VERBOSE(" LR: Generate C22 Covariance Matrix... (Prediction Stage)");
    START_TIMING(mat_gen_time);

    HICMA_problem_t hicma_problem;
    hicma_problem.theta = (double *) theta;
    hicma_problem.noise = 1e-4;
    hicma_problem.ndim = 2;

    double *xycoord = (double *) malloc(4 * nZobs * sizeof(double));
    int i;

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0) {
        int j = 0;
        for (i = 0; i < 2 * nZobs; i++) {
            xycoord[i] = data->lobs.x[j];
            xycoord[2 * nZobs + i] = data->lobs.y[j];
            if (i % 2 != 0)
                j++;
        }
        hicma_problem.kernel_type =
                strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_PARSIMONIOUS_GCD : STARSH_SPATIAL_PARSIMONIOUS_SIMD;
    } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious")) {
        int j = 0;
        for (i = 0; i < 3 * nZobs; i++) {
            xycoord[i] = data->lobs.x[j];
            xycoord[3 * nZobs + i] = data->lobs.y[j];
            if (i % 3 != 0)
                j++;
        }
        printf("STARSH_SPATIAL_TRIPARSIMONIOUS_GCD is not implemented yet\n");
        exit(0);

    } else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
               strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        int j = 0;
        for (i = 0; i < nZobs; i++) {
            if (i == nZobs / 2)
                j = 0;
            xycoord[i] = data->lobs.x[i];
            xycoord[nZobs + i] = data->lobs.y[i];
            j++;
        }
        hicma_problem.kernel_type =
                strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_PARSIMONIOUS2_GCD : STARSH_SPATIAL_PARSIMONIOUS2_SIMD;
    } else {
        // I need to change the location struct to vector array (TO DO)
        for (i = 0; i < nZobs; i++) {
            xycoord[i] = data->lobs.x[i];
            xycoord[nZobs + i] = data->lobs.y[i];

        }
        hicma_problem.kernel_type =
                strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_MATERN2_GCD : STARSH_SPATIAL_MATERN2_SIMD;
    }
    hicma_problem.point = xycoord;
    HICMA_zgenerate_problem(hicma_data_type, 'S', 0, nZobs, lts, hicma_descC22UV->mt, hicma_descC22UV->nt,
                            &hicma_problem);
    mpiF = hicma_problem.starsh_format;
    HICMA_zgytlr_Tile(ChamLower, hicma_descC22UV, hicma_descC22D, hicma_descC22rk, 0, maxrank, pow(10, -1.0 * acc),
                      compress_diag, hicma_descC);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix
    VERBOSE("LR: Generate C12 Covariance Matrix... (Prediction Stage)");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC12, &data->lmiss, &data->lobs, &data->lm, theta, data->dm,
                                   data->kernel_fun, hsequence, hrequest);
    CHAMELEON_Sequence_Wait(hsequence);

    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);

    //***************************************
    START_TIMING(time_solve);
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
    HICMA_dpotrf_Tile(ChamLower, hicma_descC22UV, hicma_descC22D, hicma_descC22rk, 0, maxrank, pow(10, -1.0 * acc));
    HICMA_dtrsmd_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, hicma_descC22UV, hicma_descC22D,
                      hicma_descC22rk, CHAMELEON_descZobs, maxrank);
    HICMA_dtrsmd_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, hicma_descC22UV, hicma_descC22D, hicma_descC22rk,
                      CHAMELEON_descZobs, maxrank);
    flops = flops + FLOPS_DPOTRF(nZobs);
    flops = flops + FLOPS_DTRSM(ChamLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);

    //**********************************dgemm
    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
    CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAMELEON_descZobs, 0, CHAMELEON_descZmiss);
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //return back descZmiss to zmiss vector
    CHAMELEON_Tile_to_Lapack(CHAMELEON_descZmiss, Zmiss, nZmiss);

    //Estimate Mean Square Error
    if (Zactual != NULL) {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
        EXAGEOSTAT_MLE_dmse_Tile_Async(CHAMELEON_descZactual, CHAMELEON_descZmiss, CHAM_descmse, hsequence, hrequest);
        VERBOSE(" Done.\n");

        CHAMELEON_Sequence_Wait(hsequence);
        STOP_TIMING(time_mse);
        data->mserror /= nZmiss;
    } else
        data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    if (data->log == 1)
        fprintf(data->pFileLog,
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, Flops: %.8f, Mean Square Error (MSE): %.8f\n\n",
                nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

    write_prediction_result("predict_result.dat", n, data->hicma_acc, 0, 0, 0, data->mserror,
                            (mat_gen_time + time_solve + time_gemm), (flops / 1e9 / (time_solve)));

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;
}


//init Hicma decriptors
void
EXAGEOSTAT_TLR_dmle_Call(MLE_data *data, int ncores, int gpus, int lts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
//! //Initiate HICMA and allocate different descriptors for
/*!  HICMA
 * Returns MLE_data data with initial values and new descriptors locations.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] ncores: number of CPU workers.
 * @param[in] gpus: number of GPU workers.
 * @param[in] dts: tile size (MB).
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * @param[in] N: number of spatial locations.
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] nZmiss: number of missing values (unknown observations).
 * */
{

    HICMA_sequence_t *msequence;
    HICMA_request_t mrequest[2] = {HICMA_SUCCESS, HICMA_SUCCESS};
    HICMA_desc_t *hicma_descC = NULL;
    HICMA_desc_t *hicma_descZ = NULL;
    CHAM_desc_t *descZ = NULL;

    CHAM_desc_t *cham_descZcpy = NULL;
    HICMA_desc_t *hicma_descZcpy = NULL;
    //CHAM_desc_t *descproduct = NULL;
    HICMA_desc_t *hicma_descdet = NULL;
    CHAM_desc_t *CHAMELEON_descZmiss = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAMELEON_descZactual = NULL;
    CHAM_desc_t *CHAMELEON_descZobs = NULL;
    HICMA_desc_t *hicma_descCD = NULL;
    HICMA_desc_t *hicma_descCUV = NULL;
    HICMA_desc_t *hicma_descCrk = NULL;
    HICMA_desc_t *hicma_descC12D = NULL;
    HICMA_desc_t *hicma_descC12UV = NULL;
    HICMA_desc_t *hicma_descC12rk = NULL;
    HICMA_desc_t *hicma_descC22D = NULL;
    HICMA_desc_t *hicma_descC22UV = NULL;
    HICMA_desc_t *hicma_descC22rk = NULL;
    int MBC, NBC, MC, NC;
    int MBD, NBD, MD, ND;
    int MBUV, NBUV, MUV, NUV;
    int MBrk, NBrk, Mrk, Nrk;

    //For ditributed system and should be removed
    double *Zcpy = (double *) malloc(N * sizeof(double));

    //CDense Descriptor
    if (data->check == 1) {
        MBC = lts;
        NBC = lts;
        MC = N;
        NC = N;
    } else {
        MBC = 1;
        NBC = 1;
        MC = lts;
        NC = lts;
    }
    HICMA_Desc_Create(&hicma_descC, NULL, HicmaRealDouble, MBC, NBC, MBC * NBC, MC, NC, 0, 0, MC, NC,
                      p_grid, q_grid);

    //CAD Descriptor
    MBD = lts;
    NBD = lts;
    MD = N;
    ND = MBD;
    HICMA_Desc_Create(&hicma_descCD, NULL, HicmaRealDouble, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND,
                      p_grid, q_grid);
    //CAD Descriptor    
    MBUV = lts;
    NBUV = 2 * data->hicma_maxrank;
    MUV = -1;
    int N_over_lts_times_lts = N / lts * lts;
    if (N_over_lts_times_lts < N) {
        MUV = N_over_lts_times_lts + lts;
    } else if (N_over_lts_times_lts == N) {
        MUV = N_over_lts_times_lts;
    } else {
        printf("%s %d: This case should not happen\n", __FILE__, __LINE__);
        exit(-1);
    }


    double expr = (double) MUV / (double) lts;
    NUV = 2 * expr * data->hicma_maxrank;

    HICMA_Desc_Create(&hicma_descCUV, NULL, HicmaRealDouble, MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0, MUV,
                      NUV, p_grid, q_grid);
    //CUV Descriptor
    MBrk = 1;
    NBrk = 1;
    Mrk = hicma_descCUV->mt;
    Nrk = hicma_descCUV->mt;

    HICMA_Desc_Create(&hicma_descCrk, NULL, HicmaRealDouble, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk,
                      Nrk, p_grid, q_grid);
    HICMA_Sequence_Create(&msequence);
    HICMA_Desc_Create(&hicma_descZ, NULL, HicmaRealDouble, lts, lts, lts * lts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&descZ, NULL, ChamRealDouble, lts, lts, lts * lts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    HICMA_Desc_Create(&hicma_descZcpy, Zcpy, HicmaRealDouble, lts, lts, lts * lts, N, 1, 0, 0, N, 1,
                      p_grid, q_grid);
    HICMA_Desc_Create(&hicma_descdet, &data->det, HicmaRealDouble, lts, lts, lts * lts, 1, 1, 0, 0, 1, 1,
                      p_grid, q_grid);

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&cham_descZcpy, Zcpy, ChamRealDouble, lts, lts, lts * lts, N, 1, 0, 0, N, 1,
                                    p_grid, q_grid);
    if (nZmiss != 0) {
        if (strcmp(data->actualZFPath, "") == 0) {
            //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZobs, &CHAMELEON_descZcpy->mat[sizeof(double)*nZmiss], ChamRealDouble, ts, ts, ts * ts, nZobs, 1, 0, 0, nZobs, 1,p_grid,q_grid);
            //CHAMELEON_descZactual=chameleon_desc_submatrix(CHAMELEON_descZcpy, 0, 0, nZmiss, 1);
            EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZobs, &Zcpy[nZmiss], ChamRealDouble, lts, lts, lts * lts,
                                            nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
            EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZactual, Zcpy, ChamRealDouble, lts, lts, lts * lts, nZmiss,
                                            1, 0, 0, nZmiss, 1, p_grid, q_grid);
        } else {
            CHAMELEON_descZobs = cham_descZcpy;
            EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZactual, NULL, ChamRealDouble, lts, lts, lts * lts, nZmiss,
                                            1, 0, 0, nZmiss, 1, p_grid, q_grid);
        }


        //C12AD Descriptor    
        MBD = lts;
        NBD = lts;
        MD = nZmiss;
        ND = MBD;
        HICMA_Desc_Create(&hicma_descC12D, NULL, HicmaRealDouble, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD,
                          ND, p_grid, q_grid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * data->hicma_maxrank;
        MUV = nZmiss;
        NUV = 2 * MUV / lts * data->hicma_maxrank;
        HICMA_Desc_Create(&hicma_descC12UV, NULL, HicmaRealDouble, MBUV, NBUV, MBUV * NBUV, MBUV, NBUV, 0,
                          0, MBUV, NBUV, p_grid, q_grid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        Mrk = hicma_descC12UV->mt;
        Nrk = hicma_descC12UV->mt;
        HICMA_Desc_Create(&hicma_descC12rk, NULL, HicmaRealDouble, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,
                          Mrk, Nrk, p_grid, q_grid);

        //C11AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZobs;
        ND = MBD;
        HICMA_Desc_Create(&hicma_descC22D, NULL, HicmaRealDouble, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD,
                          ND, p_grid, q_grid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * data->hicma_maxrank;
        MUV = nZobs;
        NUV = 2 * MUV / lts * data->hicma_maxrank;
        HICMA_Desc_Create(&hicma_descC22UV, NULL, HicmaRealDouble, MBUV, NBUV, MBUV * NBUV, MBUV, NBUV, 0,
                          0, MBUV, NBUV, p_grid, q_grid);
        //C12Ark Descriptor            
        MBrk = 1;
        NBrk = 1;
        Mrk = hicma_descC22UV->mt;
        Nrk = hicma_descC22UV->mt;
        HICMA_Desc_Create(&hicma_descC22rk, NULL, HicmaRealDouble, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,
                          Mrk, Nrk, p_grid, q_grid);


        //Other descriptors
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZmiss, NULL, ChamRealDouble, lts, lts, lts * lts, nZmiss, 1, 0,
                                        0, nZmiss, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealDouble, lts, lts, lts * lts, 1, 1, 0, 0,
                                        1, 1, p_grid, q_grid);
    }


    //Fill data struct
    data->hicma_descC = hicma_descC;
    data->hicma_descCD = hicma_descCD;
    data->hicma_descCUV = hicma_descCUV;
    data->hicma_descCrk = hicma_descCrk;
    data->hicma_descZ = hicma_descZ;
    data->hicma_descZcpy = hicma_descZcpy;
    data->hicma_descdet = hicma_descdet;
    //data->descproduct = descproduct;
    data->descZmiss = CHAMELEON_descZmiss;
    data->hicma_descC12D = hicma_descC12D;
    data->hicma_descC12UV = hicma_descC12UV;
    data->hicma_descC12rk = hicma_descC12rk;
    data->hicma_descC22D = hicma_descC22D;
    data->hicma_descC22UV = hicma_descC22UV;
    data->hicma_descC22rk = hicma_descC22rk;
    data->descmse = CHAM_descmse;
    data->descZactual = CHAMELEON_descZactual;
    data->descZobs = CHAMELEON_descZobs;
    data->descZ = descZ;
    data->hsequence = msequence;
    data->hrequest = mrequest;
    data->mserror = 0;
    //stop gsl error handler
    gsl_set_error_handler_off();
}
