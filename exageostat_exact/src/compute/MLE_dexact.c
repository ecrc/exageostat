/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_exact.c
 *
 * ExaGeoStat exact computation main functions (i.e., generate synthetic dataset, evaluate ML function, and predication).
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-23
 *
 **/
#include "../include/MLE_exact.h"

//***************************************************************************************
output results;

void EXAGEOSTAT_MLE_dzvg_Tile(MLE_data *data, double*Nrand, double*initial_theta, int n, int dts, int log)
//! Generate Observations Vector (Z) for testing Maximum
/*! Likelihood function -- CHAMELEON-sync 
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 * 	                     that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] dts: tile size (MB) is used only in the case of HiCMA not CHAMELEON.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    double*univariate_theta;
    double*univariate2_theta;
    double*univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {

        univariate_theta = (double*) malloc(3 * sizeof(double));
        univariate2_theta = (double*) malloc(3 * sizeof(double));
        univariate3_theta = (double*) malloc(3 * sizeof(double));
        univariate_theta[0] = initial_theta[0];
        univariate_theta[1] = initial_theta[2];
        univariate_theta[2] = initial_theta[3];

        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, data->descsubC11, &data->l1,
                                      &data->l1, &data->lm, univariate_theta,
                                      data->dm, "univariate_matern_stationary", msequence, &mrequest[0]);


        nu12 = 0.5 * (initial_theta[3] + initial_theta[4]);
        rho = initial_theta[5] * sqrt((tgamma(initial_theta[3] + 1) * tgamma(initial_theta[4] + 1)) /
                                      (tgamma(initial_theta[3]) * tgamma(initial_theta[4]))) *
              tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(initial_theta[0] * initial_theta[1]);

        univariate2_theta[0] = sigma_square12;
        univariate2_theta[1] = initial_theta[2];
        univariate2_theta[2] = nu12;

        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, data->descsubC12, &data->l1,
                                      &data->l1, &data->lm, univariate2_theta,
                                      data->dm, "univariate_matern_stationary", msequence, &mrequest[0]);

        univariate3_theta[0] = initial_theta[1];
        univariate3_theta[1] = initial_theta[2];
        univariate3_theta[2] = initial_theta[4];
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, data->descsubC22, &data->l1, &data->l1,
                                      &data->lm, univariate3_theta, data->dm,
                                      "univariate_matern_stationary", msequence, &mrequest[0]);
    } else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0) {
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1,
                                      &data->l1, &data->lm, initial_theta,
                                      data->dm, "univariate_matern_stationary", msequence, mrequest);
    } else {
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1,
                                      &data->l1, &data->lm, initial_theta,
                                      data->dm, data->kernel_fun, msequence, mrequest);
    }

    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
    int success = CHAMELEON_dpotrf_Tile(ChamLower, data->descC);
    //printf(" success=%d \n", success);
    //exit(0);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication    
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, data->descC, data->descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if (log == 1) {
        double*z;
        CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) (data->descZ);
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....");
#if defined(CHAMELEON_USE_MPI)
        z = (double*) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAM_descZ, z, n);
        if ( CHAMELEON_My_Mpi_Rank() == 0 )
            write_vectors(z, data, n);
        free(z);
#else
        z = CHAM_descZ->mat;
        write_vectors(z, data, n);
        //free(z);
#endif
        VERBOSE(" Done.\n");
    }

    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, data->descC);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)\n");
    VERBOSE("************************************************************\n");
}


void EXAGEOSTAT_MLE_dzcpy(MLE_data *data, double*streamdata)
//! Copy measurements vector from Lapack
/*! format to Chameleon format.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] streamdata: measurments vector in lapack format.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    VERBOSE("Copy Z from vector to decriptor.\n");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE("Done Z copying step.\n");
    VERBOSE("************************************************************\n");
}

void EXAGEOSTAT_MLE_dzvg_Tile_Async(MLE_data *data, double*Nrand, double*initial_theta, int n, int dts, int log)
//! Generate Observations Vector (Z) for testing Maximum
/*! Likelihood function -- CHAMELEON-Async
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] dts: tile size (MB) is used only in the case of HiCMA not CHAMELEON.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    //EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1,
                                  &data->l1, &data->lm, initial_theta,
                                  data->dm, data->kernel_fun, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
    int success = CHAMELEON_dpotrf_Tile_Async(ChamLower, data->descC, msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
    CHAMELEON_dtrmm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans,
                               ChamNonUnit, 1, data->descC,
                               data->descZ, msequence, mrequest);
    VERBOSE(" Done.\n");

    //if log == 1 write vector to disk
    if (log == 1) {
        double*z;
        CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) (data->descZ);
#if defined(CHAMELEON_USE_MPI)
        z = (double*) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAM_descZ, z, n);
#else
        z = CHAM_descZ->mat;
#endif
        write_vectors(z, data, n);

#if defined(CHAMELEON_USE_MPI)
        free(z);
#endif
    }

    VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous)\n");
    VERBOSE("************************************************************\n");
}


double EXAGEOSTAT_dmle_Tile(unsigned n, const double*theta, double*grad, void *CHAM_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- CHAMELEON-sync
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library. 
     * @param[in] CHAMELEON_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0, dzcpy_time = 0.0;

    int N, NRHS, success, i, num_params;
    double flops = 0.0;
    double*univariate_theta;
    double*univariate2_theta;
    double*univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    MLE_data *data = ((MLE_data *) CHAM_data);
    data->det = 0;
    data->dotp = 0;

    CHAM_desc_t *CHAM_descC = (CHAM_desc_t *) data->descC;
    CHAM_desc_t *CHAM_descsubC11 = (CHAM_desc_t *) data->descsubC11;
    CHAM_desc_t *CHAM_descsubC12 = (CHAM_desc_t *) data->descsubC12;
    CHAM_desc_t *CHAM_descsubC22 = (CHAM_desc_t *) data->descsubC22;
    CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) data->descZ;
    CHAM_desc_t *CHAM_descZ1 = (CHAM_desc_t *) data->descZ1;
    CHAM_desc_t *CHAM_descZ2 = (CHAM_desc_t *) data->descZ2;
    CHAM_desc_t *CHAM_descZcpy = (CHAM_desc_t *) data->descZcpy;
    CHAM_desc_t *CHAM_descdet = (CHAM_desc_t *) data->descdet;
    CHAM_desc_t *CHAM_descproduct = (CHAM_desc_t *) data->descproduct;
    CHAM_desc_t *CHAM_descproduct1 = (CHAM_desc_t *) data->descproduct1;
    CHAM_desc_t *CHAM_descproduct2 = (CHAM_desc_t *) data->descproduct2;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

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
    else if (strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        num_params = 7;
    else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0 ||
               strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0) 
        num_params = 10;
    else {
        fprintf(stderr, "Choosen kernel is not exist(2)!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    N = CHAM_descC->m;
    NRHS = CHAM_descZ->n;

    START_TIMING(dzcpy_time);
    if (data->iter_count == 0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAM_descZ, CHAM_descZcpy);
    if (strcmp(data->recovery_file, "") != 0 &&
        recover(data->recovery_file, data->iter_count, theta, &loglik, num_params));
    else {
        START_TIMING(dzcpy_time);
        if (data->iter_count == 0)
            //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAM_descZ, CHAM_descZcpy);
        else {
            VERBOSE("Re-store the original Z vector...");
            CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAM_descZcpy, CHAM_descZ);
            VERBOSE(" Done.\n");
        }
        STOP_TIMING(dzcpy_time);

        //Generate new co-variance matrix C based on new theta	
        VERBOSE("Generate New Covariance Matrix...");
        START_TIMING(matrix_gen_time);
        if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {

            univariate_theta = (double*) malloc(3 * sizeof(double));
            univariate2_theta = (double*) malloc(3 * sizeof(double));
            univariate3_theta = (double*) malloc(3 * sizeof(double));
            univariate_theta[0] = theta[0];
            univariate_theta[1] = theta[2];
            univariate_theta[2] = theta[3];

            EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_descsubC11, &data->l1,
                                          &data->l1, &data->lm, univariate_theta, data->dm,
                                          "univariate_matern_stationary", msequence, &mrequest[0]);


            nu12 = 0.5 * (theta[3] + theta[4]);

            rho = theta[5] * sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) /
                                  (tgamma(theta[3]) * tgamma(theta[4]))) *
                  tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(theta[0] * theta[1]);

            univariate2_theta[0] = sigma_square12;
            univariate2_theta[1] = theta[2];
            univariate2_theta[2] = nu12;
            EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_descsubC12, &data->l1,
                                          &data->l1, &data->lm, univariate2_theta,
                                          data->dm, "univariate_matern_stationary", msequence, &mrequest[0]);

            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n");

            univariate3_theta[0] = theta[1];
            univariate3_theta[1] = theta[2];
            univariate3_theta[2] = theta[4];
            EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_descsubC22, &data->l1,
                                          &data->l1, &data->lm, univariate3_theta,
                                          data->dm, "univariate_matern_stationary", msequence, &mrequest[0]);
        } else
            EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC, &data->l1,
                                          &data->l1, &data->lm, (double*) theta, data->dm,
                                          data->kernel_fun, msequence, &mrequest[0]);

        CHAMELEON_Sequence_Wait(msequence);
        STOP_TIMING(matrix_gen_time);
        VERBOSE(" Done.\n");

        VERBOSE("Cholesky factorization of Sigma...");
        START_TIMING(time_facto);
        success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_descC);
        STOP_TIMING(time_facto);
        SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
        flops = flops + FLOPS_DPOTRF(N);
        VERBOSE(" Done.\n");

        //Calculate log(|C|) --> log(square(|L|))
        VERBOSE("Calculating the log determinant ...");
        START_TIMING(logdet_calculate);
        EXAGEOSTAT_MLE_dmdet_Tile_Async(CHAM_descC, msequence, &mrequest[0], CHAM_descdet);
        CHAMELEON_Sequence_Wait(msequence);
        // printf("det: %f\n", data->det);
        logdet = 2 * data->det;
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n");

        //Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("Solving the linear system ...\n");
        START_TIMING(time_solve);
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descZ);
        STOP_TIMING(time_solve);
        flops = flops + FLOPS_DTRSM(ChamLeft, N, NRHS);
        VERBOSE(" Done.\n");

        //Calculate MLE likelihood
        VERBOSE("Calculating the MLE likelihood function ...");

        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_descproduct);

        //***************************************

        if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {

            loglik = -(N / 2) + (N / 2) * log(N) - (N / 2) * log(data->dotp) - 0.5 * logdet -
                     (double) (N / 2.0) * log(2.0 * PI);
            CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ1, CHAM_descZ1, 0, CHAM_descproduct1);
            CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ2, CHAM_descZ2, 0, CHAM_descproduct2);
            data->variance1 = (1.0 / (N / 2)) * data->dotp1;
            data->variance2 = (1.0 / (N / 2)) * data->dotp2;

        } else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0) {

            loglik = -(N / 2) + (N / 2) * log(N) - (N / 2) * log(data->dotp) - 0.5 * logdet -
                     (double) (N / 2.0) * log(2.0 * PI);

            //to be optimized
            EXAGEOSTAT_stride_vec_Tile_Async(CHAM_descZ, CHAM_descZ1, CHAM_descZ2, msequence, &mrequest[0]);
            //*********************************
            CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ1, CHAM_descZ1, 0, CHAM_descproduct1);
            CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ2, CHAM_descZ2, 0, CHAM_descproduct2);
            data->variance1 = (1.0 / (N / 2)) * data->dotp1;
            data->variance2 = (1.0 / (N / 2)) * data->dotp2;
        } else {

            loglik = -0.5 * data->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
            data->variance = theta[0];
        }
        VERBOSE(" Done.\n");

    }

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    //MPI_Bcast(theta, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif

    fprintf(stderr, " %3d- Model Parameters (", data->iter_count + 1);

    if (data->log == 1)
        fprintf(data->pFileLog, " %3d- Model Parameters (", data->iter_count + 1);

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        fprintf(stderr, "%.8f, %.8f,", data->variance1, data->variance2);
        if (data->log == 1)
            fprintf(data->pFileLog, "%.8f, %.8f,", data->variance1, data->variance2);
        i = 2;
    } else
        i = 0;
    for (; i < num_params; i++) {
        fprintf(stderr, "%.8f", theta[i]);
        if (i < num_params - 1)
            fprintf(stderr, ",");

        if (data->log == 1)
            fprintf(data->pFileLog, "%.8f, ", theta[i]);
    }

    fprintf(stderr, ")----> LogLi: %.18f\n", loglik);
    if (data->log == 1)
        fprintf(data->pFileLog, ")----> LogLi: %.18f\n", loglik);


    fprintf(stderr, " ---- Facto Time: %6.2f\n", time_facto);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter +=/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    data->final_loglik = loglik;

    //output
    results.final_loglik = loglik;
    return loglik;
}

double EXAGEOSTAT_dmle_Tile_Async(unsigned n, const double*theta, double*grad, void *CHAM_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- CHAMELEON-Async
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library.
     * @param[in] CHAMELEON_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0, dzcpy_time = 0.0, flops = 0.0;
    int N, NRHS, success;

    MLE_data *data = ((MLE_data *) CHAM_data);
    data->det = 0;
    data->dotp = 0;

    CHAM_desc_t *CHAM_descC = (CHAM_desc_t *) data->descC;
    CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) data->descZ;
    CHAM_desc_t *CHAM_descZcpy = (CHAM_desc_t *) data->descZcpy;
    CHAM_desc_t *CHAM_descdet = (CHAM_desc_t *) data->descdet;
    CHAM_desc_t *CHAM_descproduct = (CHAM_desc_t *) data->descproduct;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    N = CHAM_descC->m;
    NRHS = CHAM_descZ->n;
    START_TIMING(dzcpy_time);
    if (data->iter_count == 0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        CHAMELEON_dlacpy_Tile_Async(ChamUpperLower, CHAM_descZ, CHAM_descZcpy, msequence, mrequest);
    else {
        VERBOSE("re-store the original Z vector...");
        CHAMELEON_dlacpy_Tile_Async(ChamUpperLower, CHAM_descZcpy, CHAM_descZ, msequence, mrequest);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(dzcpy_time);

    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix...");
    START_TIMING(matrix_gen_time);
    //EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC, msequence, mrequest, &data->l1, &data->l1,(double*) theta,  data->dm);
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC, &data->l1, &data->l1, &data->lm, (double*) theta, data->dm,
                                  data->kernel_fun, msequence, &mrequest[0]);
    CHAMELEON_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma...");
    START_TIMING(time_facto);
    success = CHAMELEON_dpotrf_Tile_Async(ChamLower, CHAM_descC, msequence, mrequest);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_DPOTRF(N);
    VERBOSE(" Done.\n");

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant ...");
    START_TIMING(logdet_calculate);
    EXAGEOSTAT_MLE_dmdet_Tile_Async(CHAM_descC, msequence, &mrequest[0], CHAM_descdet);
    CHAMELEON_Sequence_Wait(msequence);
    logdet = 2 * data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system ...\n");
    START_TIMING(time_solve);
    CHAMELEON_dtrsm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descZ, msequence,
                               mrequest);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_DTRSM(ChamLeft, N, NRHS);
    VERBOSE(" Done.\n");

    //Claculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...");
    void *ws;
    CHAMELEON_dgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_descproduct, ws, msequence,
                               mrequest);
    loglik = -0.5 * data->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
    VERBOSE(" Done.\n");

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    //Print Iteration Summary
    fprintf(stderr, "------ddotproduct: %.8f ", data->dotp);
    fprintf(stderr, "------logdet: %.8f ", logdet);
    fprintf(stderr, "------expr2: %.8f ", ((double) (N / 2) * log(2 * PI)));
    //reformat
    fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness): (%.8f, %.8f, %.8f) ----> LogLi: %.8f\n",
            data->iter_count + 1, theta[0], theta[1], theta[2], loglik);

    if (data->log == 1)
        fprintf(data->pFileLog,
                " %3d- Model Parameters (variance, range, smoothness): (%.8f, %.8f, %.8f) ----> LogLi: %.8f\n",
                data->iter_count + 1, theta[0], theta[1], theta[2], loglik);

    fprintf(stderr, " ---- Facto Time: %6.2f\n", time_facto);
    fprintf(stderr, " ---- logdet Time: %6.2f\n", logdet_calculate);
    fprintf(stderr, " ---- dtrsm Time: %6.2f\n", time_solve);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter += matrix_gen_time + time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    data->final_loglik = loglik;

    return loglik;
}

void EXAGEOSTAT_dmle_Predict_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid,
                                     int mse_flag)
//! Allocate prediction operation descriptors.
/*!  
 * Returns MLE_data data with initial values and new descriptors locations.
 * @param[in] CHAMELEON_data: MLE_data struct with different MLE inputs.
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] dts: tile size (MB).
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
 * */
{
    CHAM_desc_t *CHAM_descZmiss = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *CHAM_descC22 = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAM_descmse1 = NULL;
    CHAM_desc_t *CHAM_descmse2 = NULL;
    CHAM_desc_t *CHAM_descZactual = NULL;
    CHAM_desc_t *CHAM_descZobs = NULL;
    MLE_data *data = (MLE_data *) CHAM_data;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return;
    }

    //bi-variate case
    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_flexible") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_flexible_profile") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_flexible2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_flexible2_profile") == 0) {
        nZobs *= 2;
        nZmiss *= 2;
    }

    //Descriptors Creation
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZobs, NULL, ChamRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1,
                                    p_grid, q_grid);
    if (mse_flag == 1) {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZactual, NULL, ChamRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0,
                                        nZmiss, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0,
                                        1, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse1, &data->mserror1, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0,
                                        0, 1, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse2, &data->mserror2, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0,
                                        0, 1, 1, p_grid, q_grid);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZmiss, NULL, ChamRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss,
                                    1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC12, NULL, ChamRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0,
                                    nZmiss, nZobs, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC22, NULL, ChamRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs,
                                    nZobs, p_grid, q_grid);

    //Initiate data descriptors
    data->descZmiss = CHAM_descZmiss;
    data->descC12 = CHAM_descC12;
    data->descC22 = CHAM_descC22;
    data->descmse = CHAM_descmse;
    data->descmse1 = CHAM_descmse1;
    data->descmse2 = CHAM_descmse2;
    data->descZactual = CHAM_descZactual;
    data->descZobs = CHAM_descZobs;

}


double
EXAGEOSTAT_dmle_Predict_Tile(MLE_data *CHAM_data, double*theta, int nZmiss, int nZobs, double*Zobs, double*Zactual,
                            double*Zmiss, int n)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-sync
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
    int num_params = 0;

    CHAM_desc_t *CHAM_descZmiss = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *CHAM_descC22 = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAM_descmse1 = NULL;
    CHAM_desc_t *CHAM_descmse2 = NULL;
    CHAM_desc_t *CHAM_descZactual = NULL;
    CHAM_desc_t *CHAM_descZobs = NULL;
    MLE_data *data = (MLE_data *) CHAM_data;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    data->mserror = 0;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return -1;
    }

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        data->kernel_fun = "bivariate_matern_parsimonious";
    }

    //Initiate data descriptors
    CHAM_descZmiss = data->descZmiss;
    CHAM_descC12 = data->descC12;
    CHAM_descC22 = data->descC22;
    CHAM_descmse = data->descmse;
    CHAM_descmse1 = data->descmse1;
    CHAM_descmse2 = data->descmse2;
    CHAM_descZactual = data->descZactual;
    CHAM_descZobs = data->descZobs;

    //Copy data to vectors 
    VERBOSE("Copy measurments vector to descZobs descriptor...");
    CHAMELEON_Lapack_to_Tile(Zobs, nZobs, CHAM_descZobs);
    VERBOSE(" Done.\n");

    if (Zactual != NULL) {
        VERBOSE("Copy actual measurments vector to descZactual descriptor...");
        CHAMELEON_Lapack_to_Tile(Zactual, nZmiss, CHAM_descZactual);
        VERBOSE(" Done.\n");
    }

    CHAMELEON_Sequence_Wait(msequence);

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
    else if (strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        num_params = 7;
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    printf("estimated parameters:");
    int i = 0;
    for (i = 0; i < num_params; i++) {
        printf("%.8f,", theta[i]);

    }
    printf(")\n");
    START_TIMING(mat_gen_time);

    //Generate C22 covariance matrix
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC12, &data->lmiss, &data->lobs, &data->lm, theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
    CHAMELEON_dposv_Tile(ChamLower, CHAM_descC22, CHAM_descZobs);
    flops = flops + FLOPS_DPOTRF(nZobs);
    flops = flops + FLOPS_DTRSM(ChamLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);

    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
    CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAM_descZobs, 0, CHAM_descZmiss);
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //return back descZmiss to zmiss vector
    CHAMELEON_Tile_to_Lapack(CHAM_descZmiss, Zmiss, nZmiss);

    //Estimate Mean Square Error
    if (Zactual != NULL) {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");

        if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0)

            EXAGEOSTAT_MLE_dmse_bivariate_Tile_Async(CHAM_descZactual, CHAM_descZmiss, CHAM_descmse1, CHAM_descmse2,
                                                    CHAM_descmse, msequence, mrequest);
        else

            EXAGEOSTAT_MLE_dmse_Tile_Async(CHAM_descZactual, CHAM_descZmiss, CHAM_descmse, msequence, mrequest);
        CHAMELEON_Sequence_Wait(msequence);
        VERBOSE(" Done.\n");
        STOP_TIMING(time_mse);
        data->mserror /= nZmiss;
        data->mserror1 /= (nZmiss / 2);
        data->mserror2 /= (nZmiss / 2);
    } else
        data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif

    double*z;
    if (data->log == 1)
        fprintf(data->pFileLog,
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, Flops: %.8f, Mean Square Error (MSE): %.8f\n\n",
                nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

    //output
    results.mse_pred1 = data->mserror1;
    results.mse_pred2 = data->mserror2;
    results.mse_pred = data->mserror;
    results.total_pred_time = mat_gen_time + time_solve + time_gemm;
    results.total_pred_flops = flops / 1e9 / (time_solve);

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;
}

double EXAGEOSTAT_dmle_Predict_Tile_Async(MLE_data *CHAM_data, double*theta, int nZmiss, int nZobs, int n)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-Async
 * Returns the prediction Mean Square Error (MSE) as double
 * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] n: number of spatial locations.
 * */
{

    //initialization	
    location *l1 = NULL, *l2 = NULL;
    location temp_loc;
    double mat_gen_time = 0.0;
    double time_solve = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;
    MLE_data *data = (MLE_data *) CHAM_data;
    CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) (data->descZcpy);
    CHAM_desc_t *CHAM_descZobs = (CHAM_desc_t *) (data->descZobs);
    CHAM_desc_t *CHAM_descZactual = (CHAM_desc_t *) (data->descZactual);
    CHAM_desc_t *CHAM_descZmiss = (CHAM_desc_t *) (data->descZmiss);
    CHAM_desc_t *CHAM_descC12 = (CHAM_desc_t *) (data->descC12);
    CHAM_desc_t *CHAM_descC22 = (CHAM_desc_t *) (data->descC22);
    CHAM_desc_t *CHAM_descmse = (CHAM_desc_t *) (data->descmse);
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) (data->sequence);
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    if (strcmp(data->actualZFPath, "") == 0) {
        double*z = NULL;
#if defined(CHAMELEON_USE_MPI)
        z = (double*) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAM_descZ, z, n);
#else
        z = CHAM_descZ->mat;
#endif

#if defined(CHAMELEON_USE_MPI)
        CHAMELEON_Lapack_to_Tile( z, n, CHAM_descZ);
#endif

        l1 = &data->l1;
        temp_loc.x = &l1->x[nZmiss];
        temp_loc.y = &l1->y[nZmiss];
        l2 = &temp_loc;
    } else {
        double*streamdata = NULL;
        l1 = &data->l1;
        temp_loc.x = &l1->x[nZmiss];
        temp_loc.y = &l1->y[nZmiss];
        l2 = &temp_loc;

        VERBOSE("Reading ActualZ locations for prediction from disk .....");
        l1 = readLocsFile(data->actualZLocFPath, nZmiss);
        VERBOSE(" Done.\n");

        VERBOSE("Reading ActualZ for prediction from disk .....");
        streamdata = readObsFile(data->actualZFPath, nZmiss);
        EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAM_descZactual, streamdata, msequence, mrequest);
        CHAMELEON_Sequence_Wait(data->sequence);
        VERBOSE(" Done.\n");
    }

    START_TIMING(mat_gen_time);

    //Generate C22 covariance matrix
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC22, l2, l2, &data->lm, theta, data->dm, data->kernel_fun,
                                  msequence, mrequest);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC12, l1, l2, &data->lm, theta, data->dm, data->kernel_fun,
                                  msequence, mrequest);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
    CHAMELEON_dposv_Tile_Async(ChamLower, CHAM_descC22, CHAM_descZobs, msequence, mrequest);
    flops = flops + FLOPS_DPOTRF(nZobs);
    flops = flops + FLOPS_DTRSM(ChamLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");

    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
    void *ws;
    CHAMELEON_dgemm_Tile_Async(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAM_descZobs, 0, CHAM_descZmiss, ws,
                               msequence, mrequest);
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);

    //Estimate Mean Square Error
    START_TIMING(time_mse);
    VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
    EXAGEOSTAT_MLE_dmse_Tile_Async(CHAM_descZactual, CHAM_descZmiss, CHAM_descmse, msequence, mrequest);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_mse);

    //if you do not have actual value to compare with
    if (data->descZactual == NULL)
        return -1;

    data->mserror /= nZmiss;

#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    if (data->log == 1)
        fprintf(data->pFileLog,
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, Flops: %.8f, Mean Square Error (MSE): %.8f\n\n",
                nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;
}

void EXAGEOSTAT_dmle_mloe_mmom_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid)
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

    CHAM_desc_t *CHAM_desck_t = NULL;
    CHAM_desc_t *CHAM_desck_a = NULL;
    CHAM_desc_t *CHAM_desck_atmp = NULL;
    CHAM_desc_t *CHAM_desck_ttmp = NULL;
    CHAM_desc_t *CHAM_descK_t = NULL;
    CHAM_desc_t *CHAM_descK_ttmp = NULL;
    CHAM_desc_t *CHAM_descK_a = NULL;
    CHAM_desc_t *CHAM_descexpr1 = NULL;
    CHAM_desc_t *CHAM_descexpr2 = NULL;
    CHAM_desc_t *CHAM_descexpr3 = NULL;
    CHAM_desc_t *CHAM_descexpr4 = NULL;
    CHAM_desc_t *CHAM_descmloe = NULL;
    CHAM_desc_t *CHAM_descmmom = NULL;
    CHAM_desc_t *CHAM_descalpha = NULL;
    CHAM_desc_t *CHAM_desctruthalpha = NULL;
    CHAM_desc_t *CHAM_descestimatedalpha = NULL;
    CHAM_desc_t *CHAM_desc_mloe_mmom = NULL;
    MLE_data *data = (MLE_data *) CHAM_data;
    int p = 0;
    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return;
    }
    if (strcmp(data->kernel_fun, "univariate_matern_stationary") == 0)
        p = 1;
    else if (strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        p = 1;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0)
        p = 1;
    else if (strcmp(data->kernel_fun, "bivariate_matern_flexible") == 0)
        p = 2;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0 ||
             strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        p = 2;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
             strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
        p = 2;
    else if (strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        p = 1;
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_desck_t, NULL, ChamRealDouble, dts, dts, dts * dts, p * nZobs, p, 0, 0,
                                    p * nZobs, p, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_desck_a, NULL, ChamRealDouble, dts, dts, dts * dts, p * nZobs, p, 0, 0,
                                    p * nZobs, p, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_desck_atmp, NULL, ChamRealDouble, dts, dts, dts * dts, p * nZobs, p, 0, 0,
                                    p * nZobs, p, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_desck_ttmp, NULL, ChamRealDouble, dts, dts, dts * dts, p * nZobs, p, 0, 0,
                                    p * nZobs, p, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descexpr1, NULL, ChamRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descexpr2, NULL, ChamRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descexpr3, NULL, ChamRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descexpr4, NULL, ChamRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p,
                                    p_grid, q_grid);
    data->mloe = 0;
    data->mmom = 0;
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmloe, &data->mloe, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmmom, &data->mmom, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_desctruthalpha, NULL, ChamRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descestimatedalpha, NULL, ChamRealDouble, dts, dts, dts * dts, p, p, 0, 0, p,
                                    p, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descK_t, NULL, ChamRealDouble, dts, dts, dts * dts, p * nZobs, p * nZobs, 0,
                                    0, p * nZobs, p * nZobs, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descK_a, NULL, ChamRealDouble, dts, dts, dts * dts, p * nZobs, p * nZobs, 0,
                                    0, p * nZobs, p * nZobs, p_grid, q_grid);

    //Initiae data descriptors
    data->desck_t = CHAM_desck_t;
    data->desck_a = CHAM_desck_a;
    data->desck_atmp = CHAM_desck_atmp;
    data->desck_ttmp = CHAM_desck_ttmp;
    data->descK_t = CHAM_descK_t;
    data->descK_a = CHAM_descK_a;
    data->descexpr1 = CHAM_descexpr1;
    data->descexpr2 = CHAM_descexpr2;
    data->descexpr3 = CHAM_descexpr3;
    data->descexpr4 = CHAM_descexpr4;
    data->descmloe = CHAM_descmloe;
    data->descmmom = CHAM_descmmom;
    data->descestimatedalpha = CHAM_descestimatedalpha;
    data->desctruthalpha = CHAM_desctruthalpha;
}


void
EXAGEOSTAT_dmle_mloe_mmom_Tile(MLE_data *CHAM_data, double*truth_theta, double*estimated_theta, int nZmiss, int nZobs,
                              int n)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-sync
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
    printf("%f, %f, %f, %f\n", truth_theta[0], truth_theta[1], truth_theta[2], truth_theta[3]);
    printf("%f, %f, %f,%f \n", estimated_theta[0], estimated_theta[1], estimated_theta[2], estimated_theta[3]);
    double loe_sum = 0.0;
    double mom_sum = 0.0;
    int i = 0;
    int p = 0;
    int v = 0;
    int j = 0;
    double all_time = 0.0;
    double cholesky1 = 0.0;
    double cholesky2 = 0.0;
    double matrix_gen = 0.0;
    double vecs_gen = 0.0;
    double matrix_gen3 = 0.0;
    double matrix_gen4 = 0.0;
    double copy_vecs = 0.0;
    double trsm1 = 0.0;
    double trsm2 = 0.0;
    double trsm3 = 0.0;
    double trsm4 = 0.0;
    double trsm5 = 0.0;
    double trsm6 = 0.0;
    double trsm7 = 0.0;
    double gevv1 = 0.0;
    double gevv2 = 0.0;
    double gevv3 = 0.0;
    double gevv4 = 0.0;
    double gevv5 = 0.0;
    //************************************************************************
    double*loe = (double*) malloc(nZmiss * sizeof(double));
    double*mom = (double*) malloc(nZmiss * sizeof(double));

    MLE_data *data = (MLE_data *) CHAM_data;
    CHAM_desc_t *CHAM_desck_t = data->desck_t;
    CHAM_desc_t *CHAM_desck_a = data->desck_a;
    CHAM_desc_t *CHAM_descK_t = data->descK_t;
    CHAM_desc_t *CHAM_descK_a = data->descK_a;
    CHAM_desc_t *CHAM_desck_atmp = data->desck_atmp;
    CHAM_desc_t *CHAM_desck_ttmp = data->desck_ttmp;
    CHAM_desc_t *CHAM_descexpr1 = data->descexpr1;
    CHAM_desc_t *CHAM_descexpr2 = data->descexpr2;
    CHAM_desc_t *CHAM_descexpr3 = data->descexpr3;
    CHAM_desc_t *CHAM_descexpr4 = data->descexpr4;
    CHAM_desc_t *CHAM_descmloe = data->descmloe;
    CHAM_desc_t *CHAM_descmmom = data->descmmom;
    CHAM_desc_t *CHAM_descestimatedalpha = data->descestimatedalpha;
    CHAM_desc_t *CHAM_desctruthalpha = data->desctruthalpha;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) (data->sequence);
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    location lmiss;
    lmiss.x = (double*) malloc(sizeof(double));
    lmiss.y = (double*) malloc(sizeof(double));

    double*univariate_theta;
    double*univariate2_theta;
    double*univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    double flops = 0.0;
    START_TIMING(all_time);

    int m = CHAM_descestimatedalpha->m;
    double*truthalpha = (double*) malloc(m * m * sizeof(double));
    double*estimatedalpha = (double*) malloc(m * m * sizeof(double));
    double*temp1 = (double*) malloc(m * m * sizeof(double));
    double*temp2 = (double*) malloc(m * m * sizeof(double));
    double*temp3 = (double*) malloc(m * m * sizeof(double));
    if (m == 1) {
        truthalpha[0] = truth_theta[0];
        estimatedalpha[0] = estimated_theta[0];
    }

    if (m == 2) {
        double truth_nu12 = 0.5 * (truth_theta[3] + truth_theta[4]);
        double truth_rho = truth_theta[5] * sqrt((tgamma(truth_theta[3] + 1) * tgamma(truth_theta[4] + 1)) /
                                                 (tgamma(truth_theta[3]) * tgamma(truth_theta[4]))) *
                           tgamma(truth_nu12) / tgamma(truth_nu12 + 1);

        double estimated_nu12 = 0.5 * (estimated_theta[3] + estimated_theta[4]);
        double estimated_rho =
                estimated_theta[5] * sqrt((tgamma(estimated_theta[3] + 1) * tgamma(estimated_theta[4] + 1)) /
                                          (tgamma(estimated_theta[3]) * tgamma(estimated_theta[4]))) *
                tgamma(estimated_nu12) / tgamma(estimated_nu12 + 1);

        truthalpha[0] = truth_theta[0];
        estimatedalpha[0] = estimated_theta[0];

        truthalpha[1] = truthalpha[3] = truth_rho
                                        * sqrt(truth_theta[0] * truth_theta[1]);

        estimatedalpha[1] = estimatedalpha[3] = estimated_rho
                                                * sqrt(estimated_theta[0] * estimated_theta[1]);
        truthalpha[2] = truth_theta[1];
        estimatedalpha[2] = estimated_theta[1];

    }
    CHAMELEON_Lapack_to_Tile(truthalpha, m, CHAM_desctruthalpha);
    CHAMELEON_Lapack_to_Tile(estimatedalpha, m, CHAM_descestimatedalpha);

    char* name = malloc(strlen(data->kernel_fun) + 1);
    strcpy(name, data->kernel_fun);

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        printf("TODO:  running with Z=1 no more\n");
        data->kernel_fun = "bivariate_matern_parsimonious";
    }

    START_TIMING(matrix_gen);
    VERBOSE("Create K_a and K_t Covariance Matrices (MLOE-MMOM).....");

    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descK_a, &data->lobs, &data->lobs, &data->lm, estimated_theta,
                                  data->dm, data->kernel_fun, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    //(2)Generate the co-variance matrix descK_t
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descK_t, &data->lobs, &data->lobs, &data->lm, truth_theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");
    STOP_TIMING(matrix_gen);

    //Cholesky factorization for the Co-variance matrix CHAM_descK_a
    START_TIMING(cholesky1);
    VERBOSE("(3)Cholesky factorization of CHAM_descK_a (MLOE-MMOM) .....");
    int success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_descK_a);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");
    STOP_TIMING(cholesky1);
    flops = flops + FLOPS_DPOTRF(CHAM_descK_a->m);

    START_TIMING(cholesky2);
    //(5)Cholesky factorization for the Co-variance matrix CHAM_descK_t
    VERBOSE("(5)Cholesky factorization of CHAM_descK_t (MLOE-MMOM) .....");
    success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_descK_t);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");
    STOP_TIMING(cholesky2);
    flops = flops + FLOPS_DPOTRF(CHAM_descK_t->m);

    double total_loop_time = 0.0;
    double loop_time = 0.0;
    for (p = 0; p < nZmiss; p++) {
#if defined(CHAMELEON_USE_MPI)
        if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
        fprintf(stderr, "%d ", p);
#if defined(CHAMELEON_USE_MPI)
        }
#endif
        lmiss.x[0] = data->lmiss.x[p];
        lmiss.y[0] = data->lmiss.y[p];

        VERBOSE("Generate two vectors k_a and k_t (MLOE-MMOM).....");
        START_TIMING(vecs_gen);
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_desck_t, &data->lobs, &lmiss, &data->lm, truth_theta,
                                      data->dm, data->kernel_fun, msequence, mrequest);
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_desck_a, &data->lobs, &lmiss, &data->lm, estimated_theta,
                                      data->dm, data->kernel_fun, msequence, mrequest);
        CHAMELEON_Sequence_Wait(msequence);
        STOP_TIMING(vecs_gen);

        //(6a)Copy CHAM_desck_a to CHAM_descK_atmp  (MLOE-MMOM)
        VERBOSE("(6a)Copy CHAM_desck_a to CHAM_descK_atmp  (MLOE-MMOM).....");
        START_TIMING(copy_vecs);
        CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAM_desck_t, CHAM_desck_ttmp);
        CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAM_desck_a, CHAM_desck_atmp);
        STOP_TIMING(copy_vecs);
        VERBOSE(" Done.\n");

        START_TIMING(loop_time);
        START_TIMING(trsm1);
        //(7) Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
        VERBOSE("Solving the linear system k_a = TRSM(l_a^-1, k_a) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descK_a, CHAM_desck_a);
        VERBOSE(" Done.\n");
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_a->m, CHAM_desck_a->n);
        STOP_TIMING(trsm1);

        START_TIMING(trsm2);
        //(9) Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
        VERBOSE("(9)Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_t);
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_t->m, CHAM_desck_t->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(trsm2);

        START_TIMING(trsm3);
        //(8) Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
        VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_a, CHAM_desck_a);
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_a->m, CHAM_desck_a->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(trsm3);

        START_TIMING(trsm4);
        //(10) Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
        VERBOSE("(10)Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_t);
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_t->m, CHAM_desck_t->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(trsm4);

        START_TIMING(gevv2);
        //(12) Calculate dgemm value= CHAM_desck_t^T * CHAM_desck_a
        VERBOSE("(12)Calculate dgemm CHAM_descexpr1 = CHAM_desck_t^T * CHAM_desck_a... (MLOE-MMOM)");
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desck_ttmp, CHAM_desck_a, 0, CHAM_descexpr1);
        flops = flops + FLOPS_DGEMM(CHAM_desck_ttmp->m, CHAM_desck_a->n, CHAM_descexpr1->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv2);
        START_TIMING(gevv3);
        //(13) Calculate dgemm value= CHAM_desck_a^T * CHAM_desck_atmp
        VERBOSE("(13)Calculate dgemm CHAM_descexpr1 = CHAM_desck_a^T * CHAM_desck_a... (MLOE-MMOM)");
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desck_atmp, CHAM_desck_a, 0, CHAM_descexpr4);
        flops = flops + FLOPS_DGEMM(CHAM_desck_atmp->m, CHAM_desck_a->n, CHAM_descexpr4->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv3);


        START_TIMING(gevv1);
        //(11) Calculate dgemm value= CHAM_desck_a^T * CHAM_desck_t
        VERBOSE("(11)Calculate dgemm CHAM_descexpr4 = CHAM_desck_a^T * CHAM_desck_t... (Prediction Stage)");
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desck_ttmp, CHAM_desck_t, 0, CHAM_descexpr3);
        flops = flops + FLOPS_DGEMM(CHAM_desck_ttmp->m, CHAM_desck_t->n, CHAM_descexpr3->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv1);

        //(14) Calculate dgemm CHAM_desck_a= CHAM_descK_t * CHAM_desck_a (use k_t as k_a)
        START_TIMING(gevv4);
        CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_a);
        STOP_TIMING(gevv4);

        //(13) Calculate dgemm value= CHAM_desck_a^T * CHAM_desck_t
        VERBOSE("(17)Calculate dgemm CHAM_descexpr1 = CHAM_desck_a^T * CHAM_desck_a... (Prediction Stage)");
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desck_a, CHAM_desck_a, 0, CHAM_descexpr2);
        flops = flops + FLOPS_DGEMM(CHAM_desck_atmp->m, CHAM_desck_t->n, CHAM_descexpr2->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv5);

        STOP_TIMING(loop_time);

        total_loop_time += loop_time;

        CHAMELEON_dgeadd_Tile(ChamNoTrans, 1, CHAM_desctruthalpha, -2, CHAM_descexpr1);
        CHAMELEON_dgeadd_Tile(ChamNoTrans, 1, CHAM_descexpr1, 1, CHAM_descexpr2);

        CHAMELEON_dgeadd_Tile(ChamNoTrans, 1, CHAM_desctruthalpha, -1, CHAM_descexpr3);
        CHAMELEON_dgeadd_Tile(ChamNoTrans, 1, CHAM_descestimatedalpha, -1, CHAM_descexpr4);

        printf("%f, %f, %f, %f\n", matrix_gen, vecs_gen, cholesky1, cholesky2);
        printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", trsm1, trsm2, trsm3, trsm5, trsm6, gevv1, gevv2, gevv3,
               gevv4, gevv5);
        //*******************************

        EXAGEOSTAT_MLE_dmloe_mmom_Tile_Async(CHAM_descexpr2, CHAM_descexpr3, CHAM_descexpr4, CHAM_descmloe,
                                            CHAM_descmmom, msequence, mrequest);
        CHAMELEON_Sequence_Wait(msequence);

    }
#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif

    printf("\nnZmiss = %d\n", nZmiss);
    printf(" ----(MLOE-MMOM) Gflop/s: %6.2f\n", flops / 1e9 / (total_loop_time + cholesky1 + cholesky2));


#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->mloe /= nZmiss;
    data->mmom /= nZmiss;
    STOP_TIMING(all_time);
    free(loe);
    free(mom);
    free(lmiss.x);
    free(lmiss.y);
    free(temp1);
    free(temp2);
    free(temp3);
    free(estimatedalpha);
    free(truthalpha);

    FILE *pFile;
#if defined(CHAMELEON_USE_MPI)
    if ( CHAMELEON_My_Mpi_Rank() == 0 )
    {
#endif
    //output
    results.mloe = data->mloe;
    results.mmom = data->mmom;
    printf("\nMMLOE = %f", data->mloe);
    printf("\nMMOM = %f\n\n\n", data->mmom);
    results.mloe_exec = "sync";
    results.total_mloe_mmom_time = all_time;
    results.matrix_gen_mloe_mmom_time = matrix_gen;
    results.cho_fact_mloe_mmom_time = cholesky1 + cholesky2;
    results.loop_mloe_mmom_time = total_loop_time;
    results.total_mloe_mmom_flops = flops / 1e9 / (total_loop_time + cholesky1 + cholesky2);
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    fprintf(stderr, " ---- mloe_mmom Time: %6.2f seconds\n\n", all_time);
    fprintf(stderr, " ---- gevv Time: %6.2f seconds\n\n", gevv4);
}

void EXAGEOSTAT_dmle_mloe_mmom_Tile_Async(MLE_data *CHAM_data, double*truth_theta, double*estimated_theta, int nZmiss,
                                         int nZobs, int n)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-sync
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

    printf("%f, %f, %f\n", truth_theta[0], truth_theta[1], truth_theta[2]);
    printf("%f, %f, %f\n", estimated_theta[0], estimated_theta[1], estimated_theta[2]);
    double loe_sum = 0.0;
    double mom_sum = 0.0;
    int i = 0;
    int p = 0;
    int v = 0;
    int j = 0;
    double all_time = 0.0;
    double cholesky1 = 0.0;
    double cholesky2 = 0.0;
    double matrix_gen1 = 0.0;
    double matrix_gen2 = 0.0;
    double vecs_gen = 0.0;
    double copy1 = 0.0;
    double copy2 = 0.0;
    double copy_vecs = 0.0;
    double trsm1 = 0.0;
    double trsm2 = 0.0;
    double trsm3 = 0.0;
    double trsm4 = 0.0;
    double trsm5 = 0.0;
    double trsm6 = 0.0;
    double trsm7 = 0.0;
    double gevv1 = 0.0;
    double gevv2 = 0.0;
    double gevv3 = 0.0;
    double gevv4 = 0.0;
    double gevv5 = 0.0;
    //************************************************************************
    double*loe = (double*) malloc(nZmiss * sizeof(double));
    double*mom = (double*) malloc(nZmiss * sizeof(double));

    MLE_data *data = (MLE_data *) CHAM_data;
    CHAM_desc_t *CHAM_desck_t = data->desck_t;
    CHAM_desc_t *CHAM_desck_a = data->desck_a;
    CHAM_desc_t *CHAM_descK_t = data->descK_t;
    CHAM_desc_t *CHAM_descK_a = data->descK_a;
    CHAM_desc_t *CHAM_desck_atmp = data->desck_atmp;
    CHAM_desc_t *CHAM_desck_ttmp = data->desck_ttmp;
    CHAM_desc_t *CHAM_descmloe = data->descmloe;
    CHAM_desc_t *CHAM_descmmom = data->descmmom;
    CHAM_desc_t *CHAM_descexpr1 = data->descexpr1;
    CHAM_desc_t *CHAM_descexpr2 = data->descexpr2;
    CHAM_desc_t *CHAM_descexpr3 = data->descexpr3;
    CHAM_desc_t *CHAM_descexpr4 = data->descexpr4;
    CHAM_desc_t *CHAM_descestimatedalpha = data->descestimatedalpha;
    CHAM_desc_t *CHAM_desctruthalpha = data->desctruthalpha;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) (data->sequence);
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    location lmiss;
    lmiss.x = (double*) malloc(sizeof(double));
    lmiss.y = (double*) malloc(sizeof(double));

    double*univariate_theta;
    double*univariate2_theta;
    double*univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    double flops = 0.0;
    START_TIMING(all_time);

    int m = CHAM_descestimatedalpha->m;
    double*truthalpha = (double*) malloc(m * m * sizeof(double));
    double*estimatedalpha = (double*) malloc(m * m * sizeof(double));
    double*temp1 = (double*) malloc(m * m * sizeof(double));
    double*temp2 = (double*) malloc(m * m * sizeof(double));
    double*temp3 = (double*) malloc(m * m * sizeof(double));
    if (m == 1) {
        truthalpha[0] = truth_theta[0];
        estimatedalpha[0] = estimated_theta[0];
    }

    if (m == 2) {
        double truth_nu12 = 0.5 * (truth_theta[3] + truth_theta[4]);
        double truth_rho = truth_theta[5] * sqrt((tgamma(truth_theta[3] + 1) * tgamma(truth_theta[4] + 1)) /
                                                 (tgamma(truth_theta[3]) * tgamma(truth_theta[4]))) *
                           tgamma(truth_nu12) / tgamma(truth_nu12 + 1);

        double estimated_nu12 = 0.5 * (estimated_theta[3] + estimated_theta[4]);
        double estimated_rho =
                estimated_theta[5] * sqrt((tgamma(estimated_theta[3] + 1) * tgamma(estimated_theta[4] + 1)) /
                                          (tgamma(estimated_theta[3]) * tgamma(estimated_theta[4]))) *
                tgamma(estimated_nu12) / tgamma(estimated_nu12 + 1);

        truthalpha[0] = truth_theta[0];
        estimatedalpha[0] = estimated_theta[0];

        truthalpha[1] = truthalpha[3] = truth_rho
                                        * sqrt(truth_theta[0] * truth_theta[1]);

        estimatedalpha[1] = estimatedalpha[3] = estimated_rho
                                                * sqrt(estimated_theta[0] * estimated_theta[1]);
        truthalpha[2] = truth_theta[1];
        estimatedalpha[2] = estimated_theta[1];

    }
    CHAMELEON_Lapack_to_Tile(truthalpha, m, CHAM_desctruthalpha);
    CHAMELEON_Lapack_to_Tile(estimatedalpha, m, CHAM_descestimatedalpha);


    char* name = malloc(strlen(data->kernel_fun) + 1);
    strcpy(name, data->kernel_fun);

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        printf("TODO:  running with Z=1 no more\n");
        data->kernel_fun = "bivariate_matern_parsimonious";
    }

    START_TIMING(matrix_gen1);
    VERBOSE("Create CHAM_descK_a Covariance Matrix (MLOE-MMOM).....");

    //(1)Generate the co-variance matrix descK_a
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descK_a, &data->lobs, &data->lobs, &data->lm, estimated_theta,
                                  data->dm, data->kernel_fun, msequence, mrequest);
    VERBOSE(" Done.\n");
    STOP_TIMING(matrix_gen1);

    START_TIMING(matrix_gen2);
    //(2)Generate the co-variance matrix descK_t
    VERBOSE("Create CHAM_descK_t Covariance Matrix (MLOE-MMOM).....");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descK_t, &data->lobs, &data->lobs, &data->lm, truth_theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    //CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");
    STOP_TIMING(matrix_gen2);

    START_TIMING(cholesky1);
    //(3)Cholesky factorization for the Co-variance matrix CHAM_descK_a
    VERBOSE("Cholesky factorization of CHAM_descK_a (MLOE-MMOM) .....");
    int success = CHAMELEON_dpotrf_Tile_Async(ChamLower, CHAM_descK_a, msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");
    STOP_TIMING(cholesky1);
    flops = flops + FLOPS_DPOTRF(CHAM_descK_a->m);

    START_TIMING(cholesky2);
    flops = flops + FLOPS_DPOTRF(CHAM_descK_t->m);
    //(5)Cholesky factorization for the Co-variance matrix CHAM_descK_t
    VERBOSE("Cholesky factorization of CHAM_descK_t (MLOE-MMOM) .....");
    success = CHAMELEON_dpotrf_Tile_Async(ChamLower, CHAM_descK_t, msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    STOP_TIMING(cholesky2);

    double total_loop_time = 0.0;
    double loop_time = 0.0;
    for (p = 0; p < nZmiss; p++) {

        lmiss.x[0] = data->lmiss.x[p];
        lmiss.y[0] = data->lmiss.y[p];

        VERBOSE("Generate two vectors k_a and k_t (MLOE-MMOM).....");
        START_TIMING(vecs_gen);
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_desck_t, &data->lobs, &lmiss, &data->lm, truth_theta,
                                      data->dm, data->kernel_fun, msequence, mrequest);
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_desck_a, &data->lobs, &lmiss, &data->lm, estimated_theta,
                                      data->dm, data->kernel_fun, msequence, mrequest);
        STOP_TIMING(vecs_gen);


        //(6a)Copy CHAM_desck_a to CHAM_descK_atmp  (MLOE-MMOM)
        VERBOSE("(6a)Copy CHAM_desck_a to CHAM_descK_atmp  (MLOE-MMOM).....");
        START_TIMING(copy_vecs);
        CHAMELEON_dlacpy_Tile_Async(ChamUpperLower, CHAM_desck_t, CHAM_desck_ttmp, msequence, mrequest);
        CHAMELEON_dlacpy_Tile_Async(ChamUpperLower, CHAM_desck_a, CHAM_desck_atmp, msequence, mrequest);
        STOP_TIMING(copy_vecs);
        VERBOSE(" Done.\n");

        START_TIMING(loop_time);
        START_TIMING(trsm1);
        //(7) Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
        VERBOSE("Solving the linear system k_a = TRSM(l_a^-1, k_a) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descK_a, CHAM_desck_a,
                                   msequence, mrequest);
        VERBOSE(" Done.\n");
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_a->m, CHAM_desck_a->n);
        STOP_TIMING(trsm1);

        START_TIMING(trsm2);
        //(9) Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
        VERBOSE("(9)Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_t,
                                   msequence, mrequest);
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_t->m, CHAM_desck_t->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(trsm2);

        START_TIMING(trsm3);
        //(8) Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
        VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile_Async(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_a, CHAM_desck_a,
                                   msequence, mrequest);
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_a->m, CHAM_desck_a->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(trsm3);

        START_TIMING(trsm4);
        //(10) Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
        VERBOSE("(10)Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)\n");
        CHAMELEON_dtrsm_Tile_Async(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_t,
                                   msequence, mrequest);
        flops = flops + FLOPS_DTRSM(ChamLeft, CHAM_descK_t->m, CHAM_desck_t->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(trsm4);

        START_TIMING(gevv2);
        //(12) Calculate dgemm value= CHAM_desck_t^T * CHAM_desck_a
        VERBOSE("(12)Calculate dgemm CHAM_descexpr1 = CHAM_desck_t^T * CHAM_desck_a... (MLOE-MMOM)");
        void *ws;
        CHAMELEON_dgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAM_desck_ttmp, CHAM_desck_a, 0, CHAM_descexpr1, ws,
                                   msequence, mrequest);
        flops = flops + FLOPS_DGEMM(CHAM_desck_ttmp->m, CHAM_desck_a->n, CHAM_descexpr1->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv2);
        START_TIMING(gevv3);
        //(13) Calculate dgemm value= CHAM_desck_a^T * CHAM_desck_atmp
        VERBOSE("(13)Calculate dgemm CHAM_descexpr1 = CHAM_desck_a^T * CHAM_desck_a... (MLOE-MMOM)");
        CHAMELEON_dgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAM_desck_atmp, CHAM_desck_a, 0, CHAM_descexpr4, ws,
                                   msequence, mrequest);
        flops = flops + FLOPS_DGEMM(CHAM_desck_atmp->m, CHAM_desck_a->n, CHAM_descexpr4->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv3);


        START_TIMING(gevv1);
        //(11) Calculate dgemm value= CHAM_desck_a^T * CHAM_desck_t
        VERBOSE("(11)Calculate dgemm CHAM_descexpr4 = CHAM_desck_a^T * CHAM_desck_t... (Prediction Stage)");
        CHAMELEON_dgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAM_desck_ttmp, CHAM_desck_t, 0, CHAM_descexpr3, ws,
                                   msequence, mrequest);
        flops = flops + FLOPS_DGEMM(CHAM_desck_ttmp->m, CHAM_desck_t->n, CHAM_descexpr3->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv1);

        //(14) Calculate dgemm CHAM_desck_a= CHAM_descK_t * CHAM_desck_a (use k_t as k_a)
        VERBOSE("(14)Calculate 2dtrmm CHAM_desck_a = CHAM_descK_ttmp * CHAM_desck_a... (MLOE-MMOM)");
        //CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_descK_ttmp, CHAM_desck_a, 0, CHAM_desck_t);
        START_TIMING(gevv4);
        CHAMELEON_dtrmm_Tile_Async(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_a,
                                   msequence, mrequest);
        //CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_a);
        flops = flops + FLOPS_DGEMM(CHAM_desck_ttmp->m, CHAM_desck_a->n, CHAM_desck_t->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv4);

        START_TIMING(gevv5);

        //	CHAMELEON_dtrmm_Tile_Async(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descK_t, CHAM_desck_atmp, msequence, mrequest);
        //(13) Calculate dgemm value= CHAM_desck_a^T * CHAM_desck_t
        VERBOSE("(17)Calculate dgemm CHAM_descexpr1 = CHAM_desck_a^T * CHAM_desck_a... (Prediction Stage)");
        CHAMELEON_dgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAM_desck_a, CHAM_desck_a, 0, CHAM_descexpr2, ws,
                                   msequence, mrequest);
        flops = flops + FLOPS_DGEMM(CHAM_desck_atmp->m, CHAM_desck_t->n, CHAM_descexpr2->n);
        VERBOSE(" Done.\n");
        STOP_TIMING(gevv5);


        STOP_TIMING(loop_time);

        total_loop_time += loop_time;

        CHAMELEON_dgeadd_Tile_Async(ChamNoTrans, 1, CHAM_desctruthalpha, -2, CHAM_descexpr1, msequence, mrequest);
        CHAMELEON_dgeadd_Tile_Async(ChamNoTrans, 1, CHAM_descexpr1, 1, CHAM_descexpr2, msequence, mrequest);

        CHAMELEON_dgeadd_Tile_Async(ChamNoTrans, 1, CHAM_desctruthalpha, -1, CHAM_descexpr3, msequence, mrequest);
        CHAMELEON_dgeadd_Tile_Async(ChamNoTrans, 1, CHAM_descestimatedalpha, -1, CHAM_descexpr4, msequence, mrequest);

        EXAGEOSTAT_MLE_dmloe_mmom_Tile_Async(CHAM_descexpr2, CHAM_descexpr3, CHAM_descexpr4, CHAM_descmloe,
                                            CHAM_descmmom, msequence, mrequest);
        CHAMELEON_Sequence_Wait(msequence);

    }
    //	CHAMELEON_Sequence_Wait(msequence);
#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    printf("\nnZmiss = %d\n", nZmiss);
    printf(" ----(MLOE-MMOM) Gflop/s: %6.2f\n", flops / 1e9 / (total_loop_time + cholesky1 + cholesky2));
#if defined(CHAMELEON_USE_MPI)
    }
#endif
    data->mloe /= nZmiss;
    data->mmom /= nZmiss;
    STOP_TIMING(all_time);
    free(loe);
    free(mom);
    free(lmiss.x);

    free(temp1);
    free(temp2);
    free(temp3);
    free(estimatedalpha);
    free(truthalpha);
    FILE *pFile;
#if defined(CHAMELEON_USE_MPI)
    if ( CHAMELEON_My_Mpi_Rank() == 0 )
    {
#endif
    //output
    results.mloe = data->mloe;
    results.mmom = data->mmom;
    printf("\nMMLOE = %f", data->mloe);
    printf("\nMMOM = %f\n\n\n", data->mmom);
    results.mloe_exec = "async";
    results.total_mloe_mmom_time = all_time;
    results.matrix_gen_mloe_mmom_time = matrix_gen1 + matrix_gen2;
    results.cho_fact_mloe_mmom_time = cholesky1 + cholesky2;
    results.loop_mloe_mmom_time = total_loop_time;
    results.total_mloe_mmom_flops = flops / 1e9 / (total_loop_time + cholesky1 + cholesky2);
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    fprintf(stderr, " ---- mloe_mmom Time: %6.2f seconds\n\n", all_time);
    fprintf(stderr, " ---- gevv Time: %6.2f seconds\n\n", gevv4);
}


double*EXAGEOSTAT_Fisher_Tile(MLE_data *data, int N, double*initial_theta, int dts, int p_grid, int q_grid) {
    CHAM_desc_t *CHAM_descA = NULL;
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAM_descCJ = NULL;
    CHAM_desc_t *CHAM_descresult = NULL;
    CHAM_desc_t *CHAM_descCK = NULL;
    CHAM_desc_t *CHAM_descCdiag = NULL;
    CHAM_desc_t *CHAM_descCtrace = NULL;

    double*A = (double*) calloc(9, sizeof(double));
    double time = 0.0;
    double t = 0.0;

    VERBOSE("Allocate descriptor CHAM_descC  (Fisher Matrix Generation).....");
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descA, NULL, ChamRealDouble, dts, dts, dts * dts, 3, 3, 0, 0, 3, 3, p_grid,
                                    q_grid);
    VERBOSE("...Done\n");

    VERBOSE("Allocate descriptor CHAM_descC  (Fisher Matrix Generation).....");
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC, NULL, ChamRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
                                    q_grid);
    VERBOSE("...Done\n");

    VERBOSE("Allocate descriptor CHAM_descCJ  (Fisher Matrix Generation).....");
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descCJ, NULL, ChamRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
                                    q_grid);
    VERBOSE("...Done\n");

    VERBOSE("Allocate descriptor CHAM_descCK  (Fisher Matrix Generation).....");
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descCK, NULL, ChamRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descresult, NULL, ChamRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N,
                                    p_grid, q_grid);
    VERBOSE("...Done\n");

    VERBOSE("Allocate descriptor CHAM_descCdiag  (Fisher Matrix Generation).....");
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descCdiag, NULL, ChamRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1,
                                    p_grid, q_grid);
    VERBOSE("...Done\n");

    VERBOSE("Allocate descriptor CHAM_descCtrace  (Fisher Matrix Generation).....");
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descCtrace, &t, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid,
                                    q_grid);
    VERBOSE("...Done\n");

    START_TIMING(time);
    VERBOSE("Generate covariance matrix  CHAM_descC  (Fisher Matrix Generation).....");
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC, &data->l1, &data->l1, &data->lm, initial_theta,
                                  data->dm, "univariate_matern_stationary", data->sequence, data->request);
    VERBOSE("...Done\n");
    CHAMELEON_Sequence_Wait(data->sequence);

    VERBOSE("Calculate Cholesky decomposition  (Fisher Matrix Generation).....");
    int success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_descC);
    VERBOSE("...Done\n");

    char* kernel = (double*) malloc(15 * sizeof(char));
    for (int j = 0; j < 3; j++) {
        if (j == 0) { kernel = "univariate_matern_dsigma_square"; }
        else if (j == 1) { kernel = "univariate_matern_dbeta"; }
        else if (j == 2) { kernel = "univariate_matern_dnu"; }

        VERBOSE("Generate covariance matrix  CHAM_descCJ  (Fisher Matrix Generation).....");
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_descCJ, &data->l1, &data->l1, &data->lm, initial_theta,
                                      data->dm, kernel, data->sequence, data->request);
        VERBOSE("...Done\n");
        CHAMELEON_Sequence_Wait(data->sequence);

        VERBOSE("Compute tringular solve  CHAM_descCJ  (Fisher Matrix Generation).....");
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descCJ);
        CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descCJ);
        VERBOSE("...Done\n");


        for (int k = j; k < 3; k++) {

            if (k == 0) { kernel = "univariate_matern_dsigma_square"; }

            else if (k == 1) { kernel = "univariate_matern_dbeta"; }

            else if (k == 2) { kernel = "univariate_matern_dnu"; }

            VERBOSE("Generate covariance matrix  CHAM_descCJ  (Fisher Matrix Generation).....");
            EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, CHAM_descCK, &data->l1, &data->l1, &data->lm, initial_theta,
                                          data->dm, kernel, data->sequence, data->request);
            VERBOSE("...Done\n");
            CHAMELEON_Sequence_Wait(data->sequence);

            VERBOSE("Compute tringular solve  CHAM_descCK  (Fisher Matrix Generation).....");
            CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descCK);
            CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descCK);
            VERBOSE("...Done\n");
            VERBOSE("Compute matrix-matrix multiplication  CHAM_descCK  (Fisher Matrix Generation).....");
            CHAMELEON_dgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_descCJ, CHAM_descCK, 0, CHAM_descresult);
            VERBOSE("...Done\n");


            VERBOSE("Compute the trace/diagonal of CHAM_descCK (Fisher Matrix Generation).....");
            EXAGEOSTAT_MLE_dtrace_Tile_Async(CHAM_descresult, data->sequence, data->request, CHAM_descCtrace,
                                            CHAM_descCdiag);
            VERBOSE("...Done\n");
            CHAMELEON_Sequence_Wait(data->sequence);
            A[k + 3 * j] = 0.5 * t;
            t = 0;
        }

    }

    STOP_TIMING(time);

    VERBOSE("Copy A array to descriptor CHAM_descA (Fisher Matrix Generation).....");
    CHAMELEON_Lapack_to_Tile(A, 9, CHAM_descA);
    VERBOSE("...Done\n");
    print_dmatrix("Fisher matrix", 3, 3, A, 3);

    VERBOSE("Calculate Cholesky decomposition  (Fisher Matrix Generation).....");
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', 3, A, 3);
    VERBOSE("...Done\n");


    VERBOSE("Generate Identity Matrix (I) (Fisher Matrix Generation).....");
    double*I_matrix = (double*) calloc(9, sizeof(double));
    LAPACKE_dlaset(LAPACK_COL_MAJOR, 'L', 3, 3, 0, 1, I_matrix, 3);
    VERBOSE("...Done\n");

    cblas_dtrsm(
            CblasColMajor,
            CblasLeft,
            CblasLower,
            CblasNoTrans,
            CblasNonUnit,
            3, 3, 1.0, A, 3, I_matrix, 3);

    cblas_dtrsm(
            CblasColMajor,
            CblasLeft,
            CblasLower,
            CblasTrans,
            CblasNonUnit,
            3, 3, 1.0, A, 3, I_matrix, 3);

    print_dmatrix("Fisher matrix (inverse)", 3, 3, I_matrix, 3);

    results.fisher_00 = I_matrix[0];
    results.fisher_11 = I_matrix[4];
    results.fisher_22 = I_matrix[8];
    results.fisher_time = time;

    fprintf(stderr, " ---- Time to calculate Fisher matrix: %6.4f seconds\n", time);

    double*r = calloc(10, sizeof(double));

    for (int i = 0; i < 9; i++) {
        *(r + i) = I_matrix[i];
    }

    return (r);
}


//init Chameleon descriptors
void
EXAGEOSTAT_dmle_Call(MLE_data *data, int ncores, int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
//! //Initiate CHAM and allocate different descriptors for
/*!  CHAMELEON
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

    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAM_descsubC11 = NULL;
    CHAM_desc_t *CHAM_descsubC12 = NULL;
    CHAM_desc_t *CHAM_descsubC22 = NULL;
    CHAM_desc_t *CHAM_descZ = NULL;
    CHAM_desc_t *CHAM_descZ1 = NULL;
    CHAM_desc_t *CHAM_descZ2 = NULL;
    CHAM_desc_t *CHAM_descZcpy = NULL;
    CHAM_desc_t *CHAM_descproduct = NULL;
    CHAM_desc_t *CHAM_descproduct1 = NULL;
    CHAM_desc_t *CHAM_descproduct2 = NULL;
    CHAM_desc_t *CHAM_descdet = NULL;

    // For ditributed system and should be removed
    double*Zcpy = (double*) malloc(N * sizeof(double));

    //Identifies a set of routines sharing common exception handling.
    CHAMELEON_Sequence_Create(&msequence);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC, NULL, ChamRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZ, NULL, ChamRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZcpy, Zcpy, ChamRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct, &data->dotp, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1,
                                    1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct1, &data->dotp1, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0,
                                    1, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct2, &data->dotp2, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0,
                                    1, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descdet, &data->det, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZ1, NULL, ChamRealDouble, dts, dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1,
                                    p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZ2, NULL, ChamRealDouble, dts, dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1,
                                    p_grid, q_grid);


    CHAM_descsubC11 = chameleon_desc_submatrix(CHAM_descC, 0, 0, CHAM_descC->m / 2, CHAM_descC->n / 2);
    CHAM_descsubC12 = chameleon_desc_submatrix(CHAM_descC, CHAM_descC->m / 2, 0, CHAM_descC->m / 2, CHAM_descC->n / 2);
    CHAM_descsubC22 = chameleon_desc_submatrix(CHAM_descC, CHAM_descC->m / 2, CHAM_descC->n / 2, CHAM_descC->m / 2,
                                               CHAM_descC->n / 2);

    //Fill data struct
    data->descC = CHAM_descC;
    data->descsubC11 = CHAM_descsubC11;
    data->descsubC12 = CHAM_descsubC12;
    data->descsubC22 = CHAM_descsubC22;
    data->descZ = CHAM_descZ;
    data->descZ1 = CHAM_descZ1;
    data->descZ2 = CHAM_descZ2;
    data->descZcpy = CHAM_descZcpy;
    data->descdet = CHAM_descdet;
    data->descproduct = CHAM_descproduct;
    data->descproduct1 = CHAM_descproduct1;
    data->descproduct2 = CHAM_descproduct2;
    data->sequence = msequence;
    data->request = mrequest;
    //stop gsl error handler
    gsl_set_error_handler_off();
}
