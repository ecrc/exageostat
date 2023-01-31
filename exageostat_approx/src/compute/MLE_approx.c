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
 * @date 2022-11-09
 *
 **/
#include "../include/diag.h"
#include "../include/MLE_approx.h"
//***************************************************************************************

void EXAGEOSTAT_MLE_dzvg_diag_Tile(MLE_data *data, double* Nrand, double* initial_theta, int n, int ts, int log)
//! Generate Observations Vector (Z) for testing Maximum
/*! Likelihood function -- CHAM-sync
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 * 	                     that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] ts: tile size (MB) is used only in the case of HiCMA not Chameleon.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    double* univariate_theta;
    double* univariate2_theta;
    double* univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase - diagonal approximation).....");

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {

        univariate_theta = (double* ) malloc(3 * sizeof(double));
        univariate2_theta = (double* ) malloc(3 * sizeof(double));
        univariate3_theta = (double* ) malloc(3 * sizeof(double));
        univariate_theta[0] = initial_theta[0];
        univariate_theta[1] = initial_theta[2];
        univariate_theta[2] = initial_theta[3];

        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, data->descsubC11, &data->l1, &data->l1, &data->lm,
                                      univariate_theta, data->dm, "univariate_matern_stationary", msequence,
                                      &mrequest[0]);


        nu12 = 0.5 * (initial_theta[3] + initial_theta[4]);
        rho = initial_theta[5] * sqrt((tgamma(initial_theta[3] + 1) * tgamma(initial_theta[4] + 1)) /
                                      (tgamma(initial_theta[3]) * tgamma(initial_theta[4]))) *
              tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(initial_theta[0] * initial_theta[1]);

        univariate2_theta[0] = sigma_square12;
        univariate2_theta[1] = initial_theta[2];
        univariate2_theta[2] = nu12;

        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, data->descsubC12, &data->l1, &data->l1, &data->lm,
                                      univariate2_theta, data->dm, "univariate_matern_stationary", msequence,
                                      &mrequest[0]);


        univariate3_theta[0] = initial_theta[1];
        univariate3_theta[1] = initial_theta[2];
        univariate3_theta[2] = initial_theta[4];
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamUpperLower, data->descsubC22, &data->l1, &data->l1, &data->lm,
                                      univariate3_theta, data->dm, "univariate_matern_stationary", msequence,
                                      &mrequest[0]);
    } else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0) {
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1,
                                      &data->l1, &data->lm, initial_theta,
                                      data->dm, "univariate_matern_stationary", msequence, mrequest);
    } else {
        printf("%s\n", data->kernel_fun);
        EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1,
                                      &data->l1, &data->lm, initial_theta,
                                      data->dm, data->kernel_fun, msequence, mrequest);
    }

    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase - diagonal approximation) .....");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase- diagonal approximation) .....");
    int success = CHAMELEON_dpotrf_Tile(ChamLower, data->descC);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication    
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase - diagonal approximation) .....");
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, data->descC, data->descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if (log == 1) {
        double* z;
        CHAM_desc_t *CHAMELEON_descZ = (CHAM_desc_t *) (data->descZ);
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....");
#if defined(CHAMELEON_USE_MPI)
        z = (double* ) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAMELEON_descZ, z, n);
        if ( CHAMELEON_My_Mpi_Rank() == 0 )
            write_vectors(z, data, n);
        free(z);
#else
        z = CHAMELEON_descZ->mat;
        write_vectors(z, data, n);
        free(z);
#endif
        VERBOSE(" Done.\n");
    }

    CHAMELEON_dlaset_Tile_Async(ChamUpperLower, 0, 0, data->descC, msequence, mrequest);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous- diagonal approximation)\n");
    VERBOSE("************************************************************\n");
}


void EXAGEOSTAT_MLE_dzvg_diag_Tile_Async(MLE_data *data, double* Nrand, double* initial_theta, int n, int ts, int log)
//! Generate Observations Vector (Z) for testing Maximum
/*! Likelihood function -- CHAMELEON-Async
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] ts: tile size (MB) is used only in the case of HiCMA not CHAMELEON.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    //EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, data->descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
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
    CHAMELEON_dtrmm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, data->descC, data->descZ, msequence,
                               mrequest);
    VERBOSE(" Done.\n");

    //if log == 1 write vector to disk
    if (log == 1) {
        double* z;
        CHAM_desc_t *CHAMELEON_descZ = (CHAM_desc_t *) (data->descZ);
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
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous)\n");
    VERBOSE("************************************************************\n");
}


double EXAGEOSTAT_dmle_diag_Tile(unsigned n, const double* theta, double* grad, void *CHAM_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- CHAM-sync
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library. 
     * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0, dzcpy_time = 0.0;
    int N, NRHS, success, diag_thick;
    double flops = 0.0;
    double* univariate_theta;
    double* univariate2_theta;
    double* univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;
    int num_params = 0;
    MLE_data *data = ((MLE_data *) CHAM_data);
    data->det = 0;
    data->dotp = 0;
    diag_thick = data->diag_thick;
    int i = 0;
    CHAM_desc_t *CHAM_descC = (CHAM_desc_t *) data->descC;
    CHAM_desc_t *CHAM_descsubC11 = (CHAM_desc_t *) data->descsubC11;
    CHAM_desc_t *CHAM_descsubC12 = (CHAM_desc_t *) data->descsubC12;
    CHAM_desc_t *CHAM_descsubC22 = (CHAM_desc_t *) data->descsubC22;
    CHAM_desc_t *CHAMELEON_descZ = (CHAM_desc_t *) data->descZ;
    CHAM_desc_t *CHAMELEON_descZ1 = (CHAM_desc_t *) data->descZ1;
    CHAM_desc_t *CHAMELEON_descZ2 = (CHAM_desc_t *) data->descZ2;
    CHAM_desc_t *CHAMELEON_descZ3 = (CHAM_desc_t *) data->descZ3;
    CHAM_desc_t *CHAMELEON_descZcpy = (CHAM_desc_t *) data->descZcpy;
    CHAM_desc_t *CHAM_descdet = (CHAM_desc_t *) data->descdet;
    CHAM_desc_t *CHAM_descproduct = (CHAM_desc_t *) data->descproduct;
    CHAM_desc_t *CHAM_descproduct1 = (CHAM_desc_t *) data->descproduct1;
    CHAM_desc_t *CHAM_descproduct2 = (CHAM_desc_t *) data->descproduct2;
    CHAM_desc_t *CHAM_descproduct3 = (CHAM_desc_t *) data->descproduct3;
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
    else if (strcmp(data->kernel_fun, "univariate_matern_non_gaussian") == 0 ||
             strcmp(data->kernel_fun, "univariate_exp_non_gaussian") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        num_params = 10;
    else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0
             || strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile"))
        num_params = 10;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stat") == 0)
        num_params = 8;
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    N = CHAM_descC->m;
    NRHS = CHAMELEON_descZ->n;
    START_TIMING(dzcpy_time);
    if (data->iter_count == 0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAMELEON_descZ, CHAMELEON_descZcpy);
    else {
        VERBOSE("Re-store the original Z vector...");
        CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAMELEON_descZcpy, CHAMELEON_descZ);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(dzcpy_time);

  //  CHAMELEON_dlaset_Tile_Async(ChamUpperLower, 0, 0, CHAM_descC, msequence, mrequest);
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

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma...");
    START_TIMING(time_facto);
    success = CHAM_dpotrf_diag_Tile(ChamLower, CHAM_descC, diag_thick);
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
    CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAMELEON_descZ);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_DTRSM(ChamLeft, N, NRHS);
    VERBOSE(" Done.\n");


    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...");
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ, CHAMELEON_descZ, 0, CHAM_descproduct);

    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {

	    loglik = -(N / 2) + (N / 2) * log(N) - (N / 2) * log(data->dotp) - 0.5 * logdet -
		    (double) (N / 2.0) * log(2.0 * PI);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ1, CHAMELEON_descZ1, 0, CHAM_descproduct1);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ2, CHAMELEON_descZ2, 0, CHAM_descproduct2);
	    data->variance1 = (1.0 / (N / 2)) * data->dotp1;
	    data->variance2 = (1.0 / (N / 2)) * data->dotp2;
    } else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0) {

	    loglik = -(N / 2.0) + (N / 2.0) * log(N / 2.0) - (N / 2.0) * log(data->dotp) - 0.5 * logdet -
		    (double) (N / 2.0) * log(2.0 * PI);
	    //to be optimized
	    EXAGEOSTAT_stride_vec_Tile_Async(CHAMELEON_descZ, CHAMELEON_descZ1, CHAMELEON_descZ2, msequence, &mrequest[0]);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ1, CHAMELEON_descZ1, 0, CHAM_descproduct1);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ2, CHAMELEON_descZ2, 0, CHAM_descproduct2);
	    data->variance1 = (1.0 / (N / 2.0)) * data->dotp1;
	    data->variance2 = (1.0 / (N / 2.0)) * data->dotp2;

    } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {

	    loglik = -(N / 3.0) + (N / 3.0) * log(N / 3.0) - (N / 3.0) * log(data->dotp) - 0.5 * logdet -
		    (double) (N / 3.0) * log(2.0 * PI);
	    //to be optimized
	    EXAGEOSTAT_tristride_vec_Tile_Async(CHAMELEON_descZ, CHAMELEON_descZ1, CHAMELEON_descZ2, CHAMELEON_descZ3,
			    msequence, &mrequest[0]);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ1, CHAMELEON_descZ1, 0, CHAM_descproduct1);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ2, CHAMELEON_descZ2, 0, CHAM_descproduct2);
	    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ3, CHAMELEON_descZ3, 0, CHAM_descproduct3);
	    data->variance1 = (1.0 / (N / 3.0)) * data->dotp1;
	    data->variance2 = (1.0 / (N / 3.0)) * data->dotp2;
	    data->variance3 = (1.0 / (N / 3.0)) * data->dotp3;
    } else {
	    loglik = -0.5 * data->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
	    data->variance = theta[0];
    }
    VERBOSE(" Done.\n");
    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    // MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
	    printf(" %3d- Model Parameters (", data->iter_count + 1);

	    if (data->log == 1)
		    fprintf(data->pFileLog, " %3d- Model Parameters (", data->iter_count + 1);

	    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
			    strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
		    printf("%.8f, %.8f,", data->variance1, data->variance2);
		    if (data->log == 1)
			    fprintf(data->pFileLog, "%.8f, %.8f,", data->variance1, data->variance2);
		    i = 2;
	    } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {
		    printf("%.8f, %.8f, %.8f,", data->variance1, data->variance2, data->variance3);
		    if (data->log == 1)
			    fprintf(data->pFileLog, "%.8f, %.8f, %.8f,", data->variance1, data->variance2, data->variance3);
		    i = 3;
		    results.estimated_theta[0] = data->variance1;
		    results.estimated_theta[1] = data->variance2;
		    results.estimated_theta[2] = data->variance3;
	    } else
		    i = 0;
	    for (; i < num_params; i++) {
		    printf("%.8f", theta[i]);
		    if (i < num_params - 1)
			    printf(",");

		    if (data->log == 1)
			    fprintf(data->pFileLog, "%.8f, ", theta[i]);
	    }

	    printf(")----> LogLi: %.18f\n", loglik);
	    if (data->log == 1)
		    fprintf(data->pFileLog, ")----> LogLi: %.18f\n", loglik);


	    printf(" ---- Facto Time: %6.2f\n", time_facto);
	    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
	    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
	    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
	    printf(" ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);
	    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + time_solve));

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

double EXAGEOSTAT_dmle_diag_Tile_Async(unsigned n, const double* theta, double* grad, void *CHAM_data) {
	//! Maximum Likelihood Evaluation (MLE)
	/*!  -- CHAMELEON-Async
	 * Returns the loglikelihhod value for the given theta.
	 * @param[in] n: unsigned variable used by NLOPT library.
	 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
	 *                           that is used to to generate the Covariance Matrix.
	 * @param[in] grad: double variable used by NLOPT library.
	 * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
	 * */
	//Initialization
	double loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0, dzcpy_time = 0.0, flops = 0.0;
	int N, NRHS, success, diag_thick;

	MLE_data *data = ((MLE_data *) CHAM_data);
	data->det = 0;
	data->dotp = 0;
	diag_thick = data->diag_thick;
	CHAM_desc_t *CHAM_descC = (CHAM_desc_t *) data->descC;
	CHAM_desc_t *CHAMELEON_descZ = (CHAM_desc_t *) data->descZ;
	CHAM_desc_t *CHAMELEON_descZcpy = (CHAM_desc_t *) data->descZcpy;
	CHAM_desc_t *CHAM_descdet = (CHAM_desc_t *) data->descdet;
	CHAM_desc_t *CHAM_descproduct = (CHAM_desc_t *) data->descproduct;
	RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
	RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

	N = CHAM_descC->m;
	NRHS = CHAMELEON_descZ->n;
	START_TIMING(dzcpy_time);
	if (data->iter_count == 0)
		//Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
		CHAMELEON_dlacpy_Tile_Async(ChamUpperLower, CHAMELEON_descZ, CHAMELEON_descZcpy, msequence, mrequest);
	else {
		VERBOSE("re-store the original Z vector...");
		CHAMELEON_dlacpy_Tile_Async(ChamUpperLower, CHAMELEON_descZcpy, CHAMELEON_descZ, msequence, mrequest);
		VERBOSE(" Done.\n");
	}
	STOP_TIMING(dzcpy_time);


	//Generate new co-variance matrix C based on new theta	
	VERBOSE("Generate New Covariance Matrix...");
	START_TIMING(matrix_gen_time);
	CHAM_MLE_dcmg_diag_Tile_Async(ChamLower, CHAM_descC, &data->l1, &data->l1, &data->lm, (double* ) theta, data->dm,
			data->kernel_fun, diag_thick, msequence, &mrequest[0]);


	STOP_TIMING(matrix_gen_time);
	VERBOSE(" Done.\n");
	CHAMELEON_Sequence_Wait(msequence);

	//Calculate Cholesky Factorization (C=LL-1)
	VERBOSE("Cholesky factorization of Sigma...");
	START_TIMING(time_facto);
	success = CHAM_dpotrf_diag_Tile_Async(ChamLower, CHAM_descC, diag_thick, msequence, mrequest);
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
	CHAMELEON_dtrsm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAMELEON_descZ, msequence,
			mrequest);
	STOP_TIMING(time_solve);
	flops = flops + FLOPS_DTRSM(ChamLeft, N, NRHS);
	VERBOSE(" Done.\n");

	//Claculate MLE likelihood
	VERBOSE("Calculating the MLE likelihood function ...");
	void *ws;
	CHAMELEON_dgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAMELEON_descZ, CHAMELEON_descZ, 0, CHAM_descproduct, &ws,
			msequence, mrequest);
	loglik = -0.5 * data->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
	VERBOSE(" Done.\n");

	//Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
	MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	if(CHAMELEON_My_Mpi_Rank() == 0)
	{
#endif
		fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n",
				data->iter_count + 1, theta[0], theta[1], theta[2], loglik);

		if (data->log == 1)
			fprintf(data->pFileLog,
					" %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n",
					data->iter_count + 1, theta[0], theta[1], theta[2], loglik);

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

double EXAGEOSTAT_dmle_diag_Predict_Tile(MLE_data *CHAM_data, double* theta, int nZmiss, int nZobs, double* Zobs,
		double* Zactual, double* Zmiss, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- CHAMELEON-sync
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
	double time_solve = 0.0;
	double time_gemm = 0.0;
	double mat_gen_time = 0.0;
	double time_mse = 0.0;
	double flops = 0.0;

	CHAM_desc_t *CHAMELEON_descZmiss = NULL;
	CHAM_desc_t *CHAM_descC12 = NULL;
	CHAM_desc_t *CHAM_descC22 = NULL;
	CHAM_desc_t *CHAM_descmse = NULL;
	CHAM_desc_t *CHAMELEON_descZactual = NULL;
	CHAM_desc_t *CHAMELEON_descZobs = NULL;
	MLE_data *data = (MLE_data *) CHAM_data;
	RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
	RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
	data->mserror = 0;

	if (nZmiss <= 0) {
		fprintf(stderr, " Number of missing values should be positive value\n");
		return -1;
	}

	//Initiate data descriptors
	CHAMELEON_descZmiss = data->descZmiss;
	CHAM_descC12 = data->descC12;
	CHAM_descC22 = data->descC22;
	CHAM_descmse = data->descmse;
	CHAMELEON_descZactual = data->descZactual;
	CHAMELEON_descZobs = data->descZobs;

	int diag_thick = data->diag_thick;


	CHAMELEON_dlaset_Tile_Async(ChamUpperLower, 0, 0, CHAM_descC12, msequence, mrequest);
	CHAMELEON_dlaset_Tile_Async(ChamUpperLower, 0, 0, CHAM_descC22, msequence, mrequest);

	//Copy data to vectors
	VERBOSE("Copy measurments vector to descZobs descriptor...");
	CHAMELEON_Lapack_to_Tile(Zobs, nZobs, CHAMELEON_descZobs);
	VERBOSE(" Done.\n");

	if (Zactual != NULL) {
		VERBOSE("Copy actual measurments vector to descZactual descriptor...");
		//EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAMELEON_descZactual, Zactual, msequence, mrequest);
		CHAMELEON_Lapack_to_Tile(Zactual, nZmiss, CHAMELEON_descZactual);
		VERBOSE(" Done.\n");
	}

	CHAMELEON_Sequence_Wait(msequence);

#if defined(CHAMELEON_USE_MPI)
	MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
	theta[0] = data->variance;
	printf("estimated parameters: %f - %f - %f\n", theta[0], theta[1], theta[2]);


	START_TIMING(mat_gen_time);
	//Generate C22 covariance matrix
	VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
	CHAM_MLE_dcmg_diag_Tile_Async(ChamLower, CHAM_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm,
			data->kernel_fun, data->diag_thick, msequence, mrequest);

	CHAMELEON_Sequence_Wait(msequence);
	VERBOSE(" Done.\n");

	//Generate C12 covariance matrix
	VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");

	CHAM_MLE_dcmg_diag_Tile_Async(ChamLower, CHAM_descC22, &data->lmiss, &data->lobs, &data->lm, theta, data->dm,
			data->kernel_fun, diag_thick, msequence, mrequest);
	CHAMELEON_Sequence_Wait(msequence);
	VERBOSE(" Done.\n");
	STOP_TIMING(mat_gen_time);


	START_TIMING(time_solve);
	//Start prediction
	VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
	CHAM_dpotrf_diag_Tile(ChamLower, CHAM_descC22, diag_thick);
	CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC22, CHAMELEON_descZobs);
	CHAMELEON_dtrsm_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, CHAM_descC22, CHAMELEON_descZobs);
	flops = flops + FLOPS_DPOTRF(nZobs);
	flops = flops + FLOPS_DTRSM(ChamLeft, nZobs, nZobs);
	STOP_TIMING(time_solve);


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
		EXAGEOSTAT_MLE_dmse_Tile_Async(CHAMELEON_descZactual, CHAMELEON_descZmiss, CHAM_descmse, msequence, mrequest);
		VERBOSE(" Done.\n");
		CHAMELEON_Sequence_Wait(msequence);
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
					"\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n",
					nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

		write_prediction_result("predict_result.dat", n, diag_thick, 0, 0, 0, data->mserror,
				(mat_gen_time + time_solve + time_gemm), (flops / 1e9 / (time_solve)));

#if defined(CHAMELEON_USE_MPI)
	}
#endif

	return data->mserror;
}


double CHAMELEON_dmle_Predict_diag_Tile_Async(MLE_data *CHAM_data, double* theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- CHAMELEON-Async
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
	double* streamdata = NULL;
	double mat_gen_time = 0.0;
	double time_solve = 0.0;
	double time_mse = 0.0;
	double flops = 0.0;
	MLE_data *data = (MLE_data *) CHAM_data;
	CHAM_desc_t *CHAMELEON_descZ = (CHAM_desc_t *) (data->descZcpy);
	CHAM_desc_t *CHAMELEON_descZobs = (CHAM_desc_t *) (data->descZobs);
	CHAM_desc_t *CHAMELEON_descZactual = (CHAM_desc_t *) (data->descZactual);
	CHAM_desc_t *CHAMELEON_descZmiss = (CHAM_desc_t *) (data->descZmiss);
	CHAM_desc_t *CHAM_descC12 = (CHAM_desc_t *) (data->descC12);
	CHAM_desc_t *CHAM_descC22 = (CHAM_desc_t *) (data->descC22);
	CHAM_desc_t *CHAM_descmse = (CHAM_desc_t *) (data->descmse);
	RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) (data->sequence);
	RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

	if (strcmp(data->actualZFPath, "") == 0) {
		double* z = NULL;
#if defined(CHAMELEON_USE_MPI)
		z = (double* ) malloc(n * sizeof(double));
		CHAMELEON_Tile_to_Lapack( CHAMELEON_descZ, z, n);
#else
		z = CHAMELEON_descZ->mat;
#endif

		shuffle(z, &data->l1, n);

#if defined(CHAMELEON_USE_MPI)
		CHAMELEON_Lapack_to_Tile( z, n, CHAMELEON_descZ);
#endif

		l1 = &data->l1;
		temp_loc.x = &l1->x[nZmiss];
		temp_loc.y = &l1->y[nZmiss];
		l2 = &temp_loc;
	} else {
		temp_loc.x = &l1->x[nZmiss];
		temp_loc.y = &l1->y[nZmiss];
		l2 = &temp_loc;

		VERBOSE("Reading ActualZ locations for prediction from disk .....");
		l1 = readLocsFile(data->actualZLocFPath, nZmiss);
		VERBOSE(" Done.\n");

		//streamdata=(double* ) malloc(nZmiss * sizeof(double));
		VERBOSE("Reading ActualZ for prediction from disk .....");
		streamdata = readObsFile(data->actualZFPath, nZmiss);
		EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAMELEON_descZactual, streamdata, msequence, mrequest);
		CHAMELEON_Sequence_Wait(data->sequence);
		VERBOSE(" Done.\n");
	}

	START_TIMING(mat_gen_time);

	//Generate C22 covariance matrix
	VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
	CHAM_MLE_dcmg_diag_Tile_Async(ChamLower, CHAM_descC22, l2, l2, &data->lm, theta, data->dm, data->kernel_fun,
			data->diag_thick, msequence, mrequest);
	VERBOSE(" Done.\n");

	//Generate C12 covariance matrix
	VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
	CHAM_MLE_dcmg_diag_Tile_Async(ChamLower, CHAM_descC12, l1, l2, &data->lm, theta, data->dm, data->kernel_fun,
			data->diag_thick, msequence, mrequest);

	VERBOSE(" Done.\n");
	STOP_TIMING(mat_gen_time);

	START_TIMING(time_solve);
	//Start prediction
	VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
	CHAMELEON_dposv_Tile_Async(ChamLower, CHAM_descC22, CHAMELEON_descZobs, msequence, mrequest);
	flops = flops + FLOPS_DPOTRF(nZobs);
	flops = flops + FLOPS_DTRSM(ChamLeft, nZobs, nZobs);
	VERBOSE(" Done.\n");

	VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
	void *ws;
	CHAMELEON_dgemm_Tile_Async(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAMELEON_descZobs, 0, CHAMELEON_descZmiss,
			&ws, msequence, mrequest);
	flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_solve);


	//Estimate Mean Square Error
	START_TIMING(time_mse);
	VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
	EXAGEOSTAT_MLE_dmse_Tile_Async(CHAMELEON_descZactual, CHAMELEON_descZmiss, CHAM_descmse, msequence, mrequest);
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
					"\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n",
					nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

#if defined(CHAMELEON_USE_MPI)
	}
#endif

	return data->mserror;
}

//init Chameleon descriptors
void EXAGEOSTAT_dmle_diag_Call(MLE_data *data, int ncores, int gpus, int dts, int p_grid, int q_grid, int N, int nZobs,
		int nZmiss)
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
	CHAM_desc_t *CHAMELEON_descZ = NULL;
	CHAM_desc_t *CHAMELEON_descZcpy = NULL;
	CHAM_desc_t *CHAM_descproduct = NULL;
	CHAM_desc_t *CHAM_descproduct1 = NULL;
	CHAM_desc_t *CHAM_descproduct2 = NULL;
	CHAM_desc_t *CHAM_descproduct3 = NULL;
	CHAM_desc_t *CHAM_descdet = NULL;

	// For ditributed system and should be removed
	double* Zcpy = (double* ) malloc(N * sizeof(double));

	//Identifies a set of routines sharing common exception handling.
	CHAMELEON_Sequence_Create(&msequence);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC, NULL, ChamRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
			q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZ, NULL, ChamRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1,
			p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZcpy, Zcpy, ChamRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1,
			p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct, &data->dotp, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1,
			1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct1, &data->dotp1, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0,
			1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct2, &data->dotp2, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0,
			1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct3, &data->dotp3, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0,
			1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descdet, &data->det, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
			p_grid, q_grid);

	//Fill data struct
	data->descC = CHAM_descC;
	data->descZ = CHAMELEON_descZ;
	data->descZcpy = CHAMELEON_descZcpy;
	data->descdet = CHAM_descdet;
	data->descproduct = CHAM_descproduct;
	data->descproduct1 = CHAM_descproduct;
	data->descproduct2 = CHAM_descproduct;
	data->descproduct3 = CHAM_descproduct;
	data->sequence = msequence;
	data->request = mrequest;
	//stop gsl error handler
	gsl_set_error_handler_off();
}
