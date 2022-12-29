/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#include "../include/MLE_exact.h"
//***************************************************************************************
void MORSE_MLE_dzvg_Tile (MLE_data *data,  double * Nrand, double * initial_theta, int n, int dts, int log)
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- MORSE-sync 
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     * 	                     that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] dts: tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * */
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    double* univariate_theta;
    double* univariate2_theta;
    double* univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;
    double time_facto = 0.0, time_trmm = 0.0, matrix_gen_time = 0.0;
    double flops = 0;
    MORSE_desc_t *MORSE_descZ = (MORSE_desc_t *)(data->descZ);
    //In the case of testing mode, Z should be generated using Nrand and initial_theta
    //if (test == 1)    
    //{
    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    START_TIMING(matrix_gen_time);
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
    {

        univariate_theta =(double *) malloc(3 * sizeof(double));
        univariate2_theta =(double *) malloc(3 * sizeof(double));
        univariate3_theta =(double *) malloc(3 * sizeof(double));
        univariate_theta[0]=initial_theta[0];
        univariate_theta[1]=initial_theta[2];
        univariate_theta[2]=initial_theta[3];

        MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, data->descsubC11, &data->l1, 
                &data->l1, &data->lm, univariate_theta, 
                data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);


        nu12 = 0.5 * (initial_theta[3] + initial_theta[4]);
        rho = initial_theta[5] * sqrt( (tgamma(initial_theta[3] + 1)*tgamma(initial_theta[4] + 1)) /
                (tgamma(initial_theta[3]) * tgamma(initial_theta[4])) ) *
            tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(initial_theta[0]*initial_theta[1]) ;

        univariate2_theta[0]=sigma_square12;
        univariate2_theta[1]=initial_theta[2];
        univariate2_theta[2]=nu12;

        MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, data->descsubC12, &data->l1, 
                &data->l1, &data->lm, univariate2_theta, 
                data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

        //	MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descsubC21, &data->l1, &data->l1, &data->lm, univariate2_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

        univariate3_theta[0]=initial_theta[1];
        univariate3_theta[1]=initial_theta[2];
        univariate3_theta[2]=initial_theta[4];
        MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, data->descsubC22, &data->l1, &data->l1, 
                &data->lm, univariate3_theta, data->dm,
                "univariate_matern_stationary",  msequence, &mrequest[0]);
    }
    else	if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
    {
        MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, &data->l1,
                &data->l1, &data->lm, initial_theta,
                data->dm, "univariate_matern_stationary",  msequence, mrequest);
    }

    else
    {
        printf("%s\n", data->kernel_fun);
        MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, &data->l1,
                &data->l1, &data->lm, initial_theta,
                data->dm, data->kernel_fun,  msequence, mrequest);
    }

    MORSE_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");
    //                MORSE_Sequence_Wait(msequence);
    //double sum=0;
    //    double *C = (double *) malloc(n * n * sizeof(double));
    //      MORSE_Tile_to_Lapack( data->descC, C, n);

    //int i=0;
    //for(i=0;i<n*n;i++)
    //	sum+=C[i];
    //printf("sum= %f\n", sum);
    //exit(0);
    //       print_dmatrix("testC", 16, 16, C, 16);
    //exit(0);

    //double *C = (double *) malloc(n *n *sizeof(double));
    //	MORSE_Tile_to_Lapack( data->descC, C, n);
    //	print_dmatrix("testC", 16, 16, C, 16);
    //	exit(0);


    //	 double *C = (double *) malloc(n * n * sizeof(double));
    //	 MORSE_Tile_to_Lapack( data->descC, C, n);
    //	 print_dmatrix("testC", n, n, C, n);
    //exit(0);
    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
    START_TIMING(time_facto);    
    int success = MORSE_dpotrf_Tile(MorseLower, data->descC);
    STOP_TIMING(time_facto);
    flops = flops + FLOPS_DPOTRF(n);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");
    //double *C = (double *) malloc(n * n * sizeof(double));
    //MORSE_Tile_to_Lapack( data->descC, C, n);
    //print_dmatrix("testC", 16, 16, C, 16);
    //exit(0);

    //Triangular matrix-matrix multiplication    
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
    START_TIMING(time_trmm);
    MORSE_dtrmm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ);
    STOP_TIMING(time_trmm);    
    flops = flops + FLOPS_DTRMM(MorseLeft, n, MORSE_descZ->n);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if(log==1)
    {
        double *z;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....");
        z = (double *) malloc(n * sizeof(double));
        MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
        if ( MORSE_My_Mpi_Rank() == 0 )
            write_vectors(z, data, n);
        free(z);
        VERBOSE(" Done.\n");
    }

    /*	}
        else
        {
        double * streamdata;
        streamdata=(double *) malloc(n * sizeof(double));

    //Reading Observations from disk and copy it to data->descZ
    VERBOSE("Reading Observations from disk .....");
    streamdata = readObsFile(data->obsFPath, n);        
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(data->sequence);
    VERBOSE(" Done.\n");
    free(streamdata);
    }
    */

    MORSE_dlaset_Tile(MorseUpperLower, 0, 0, data->descC);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)\n");
    VERBOSE("************************************************************\n");
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- dtrmm Time: %6.2f\n", time_trmm);
    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    printf(" ---- Total Time: %6.2f\n", /*matrix_gen_time +*/ time_facto + time_trmm);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + time_trmm));
    data->avg_exec_time_gen_stage = matrix_gen_time + time_facto + time_trmm;
    data->avg_flops_gen_stage = flops / 1e9 / (time_facto +time_trmm);
}



void MORSE_MLE_dzcpy( MLE_data *data, double *streamdata)
    //! Copy measurements vector from Lapack
    /*! format to Chameleon format.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] streamdata: measurments vector in lapack format.
     * */
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    VERBOSE("Copy Z from vector to decriptor.\n");
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    VERBOSE("Done Z copying step.\n");
    VERBOSE("************************************************************\n");
}

void MORSE_MLE_dzvg_Tile_Async(MLE_data *data,  double * Nrand, double * initial_theta, int n, int dts, int log)
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- MORSE-Async
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] dts: tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * */
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    //In the case of testing mode, Z should be generated using Nrand and initial_theta
    //       if (test ==1)
    //      {
    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, &data->l1, 
            &data->l1, &data->lm, initial_theta, 
            data->dm, data->kernel_fun,  msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
    int success = MORSE_dpotrf_Tile_Async(MorseLower, data->descC, msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
    MORSE_dtrmm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, 
            MorseNonUnit, 1, data->descC, 
            data->descZ, msequence, mrequest);
    VERBOSE(" Done.\n");

    //if log == 1 write vector to disk
    if(log == 1)
    {
        double *z;
        MORSE_desc_t *MORSE_descZ = (MORSE_desc_t *)(data->descZ);
#if defined(CHAMELEON_USE_MPI)
        z = (double *) malloc(n * sizeof(double));
        MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
#else
        z = MORSE_descZ->mat;
#endif
        write_vectors(z, data, n);

#if defined(CHAMELEON_USE_MPI)
        free(z);
#endif
    }

    /*       }
             else
             {
             double * streamdata;
             streamdata=(double *) malloc(n * sizeof(double));

    //Reading Observations from disk and copy it to data->descZ
    VERBOSE("Reading Observations from disk .....");
    streamdata = readObsFile(data->obsFPath, n);
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(data->sequence);
    VERBOSE(" Done.\n");
    free(streamdata);
    }
    */
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous)\n");
    VERBOSE("************************************************************\n");
}



double MORSE_dmle_Tile(unsigned n, const double * theta, double * grad, void * MORSE_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- MORSE-sync
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library. 
     * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0, dzcpy_time=0.0;

    int N, NRHS, success, i, num_params;
    double flops =0.0;	
    double* univariate_theta;
    double* univariate2_theta;
    double* univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;

    MLE_data* data	= ((MLE_data*)MORSE_data);
    data->det	= 0;
    data->dotp	= 0;

    MORSE_desc_t *MORSE_descC	 = (MORSE_desc_t *) data->descC;
    MORSE_desc_t *MORSE_descsubC11   = (MORSE_desc_t *) data->descsubC11;
    MORSE_desc_t *MORSE_descsubC12   = (MORSE_desc_t *) data->descsubC12;
    //MORSE_desc_t *MORSE_descsubC21   = (MORSE_desc_t *) data->descsubC21;
    MORSE_desc_t *MORSE_descsubC22   = (MORSE_desc_t *) data->descsubC22;
    MORSE_desc_t *MORSE_descZ	 = (MORSE_desc_t *) data->descZ;
    MORSE_desc_t *MORSE_descZ1       = (MORSE_desc_t *) data->descZ1;
    MORSE_desc_t *MORSE_descZ2       = (MORSE_desc_t *) data->descZ2;
    MORSE_desc_t *MORSE_descZcpy	 = (MORSE_desc_t *) data->descZcpy; 
    MORSE_desc_t *MORSE_descdet	 = (MORSE_desc_t *) data->descdet;
    MORSE_desc_t *MORSE_descproduct	 = (MORSE_desc_t *) data->descproduct;
    MORSE_desc_t *MORSE_descproduct1 = (MORSE_desc_t *) data->descproduct1;
    MORSE_desc_t *MORSE_descproduct2 = (MORSE_desc_t *) data->descproduct2;
    MORSE_sequence_t *msequence	 = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest	 = (MORSE_request_t *) data->request;

    if(strcmp(data->kernel_fun, "univariate_matern_stationary")   == 0 || strcmp(data->kernel_fun, "univariate_pow_exp_stationary")   == 0 )
        num_params = 3;
    else if(strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary")   == 0)
        num_params = 4;
    else if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
        num_params = 9;
    else if(strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0)
        num_params = 11;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary")   == 0)
        num_params = 7;
    else
    {
        fprintf(stderr,"Choosen kernel is not exist(2)!\n");
        fprintf(stderr, "Called function is: %s\n",__func__);
        exit(0);
    }
    N	= MORSE_descC->m;
    NRHS	= MORSE_descZ->n;



    START_TIMING(dzcpy_time);
    if(data->iter_count==0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descZ, MORSE_descZcpy);
    if(strcmp(data->recovery_file,"") != 0 && recover(data->recovery_file, data->iter_count, theta, &loglik, num_params));
    else
    {
        START_TIMING(dzcpy_time);
        if(data->iter_count==0)
            //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descZ, MORSE_descZcpy); 
        else
        {	
            VERBOSE("Re-store the original Z vector...");
            MORSE_dlacpy_Tile(MorseUpperLower ,MORSE_descZcpy,MORSE_descZ);
            VERBOSE(" Done.\n");
        }
        STOP_TIMING(dzcpy_time);	

        // double *C = (double *) malloc(N * N * sizeof(double));

        //Generate new co-variance matrix C based on new theta	
        VERBOSE("Generate New Covariance Matrix...");
        START_TIMING(matrix_gen_time);	
        //MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, msequence, &mrequest[0], &data->l1, &data->l1, (double *)theta,  data->dm); 
        if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        {

            univariate_theta =(double *) malloc(3 * sizeof(double));
            univariate2_theta =(double *) malloc(3 * sizeof(double));
            univariate3_theta =(double *) malloc(3 * sizeof(double));
            univariate_theta[0]=theta[0];
            univariate_theta[1]=theta[2];
            univariate_theta[2]=theta[3];

            MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descsubC11, &data->l1,
                    &data->l1, &data->lm, univariate_theta, data->dm,
                    "univariate_matern_stationary",  msequence, &mrequest[0]);


            nu12 = 0.5 * (theta[3] + theta[4]);

            rho = theta[5] * sqrt( (tgamma(theta[3] + 1)*tgamma(theta[4] + 1)) /
                    (tgamma(theta[3]) * tgamma(theta[4])) ) *
                tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(theta[0]*theta[1]) ;

            univariate2_theta[0]=sigma_square12;
            univariate2_theta[1]=theta[2];
            univariate2_theta[2]=nu12; 
            MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descsubC12, &data->l1,
                    &data->l1, &data->lm, univariate2_theta, 
                    data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

            //MORSE_Sequence_Wait(msequence);
            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n");

            univariate3_theta[0]=theta[1];
            univariate3_theta[1]=theta[2];
            univariate3_theta[2]=theta[4];
            MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descsubC22, &data->l1, 
                    &data->l1, &data->lm, univariate3_theta, 
                    data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);
            //MORSE_Sequence_Wait(msequence);
        }
        else   
            MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, &data->l1,
                    &data->l1, &data->lm, (double *)theta, data->dm,
                    data->kernel_fun,  msequence, &mrequest[0]);

        MORSE_Sequence_Wait(msequence);
        STOP_TIMING(matrix_gen_time);
        VERBOSE(" Done.\n");

        // double *C = (double *) malloc(N * N * sizeof(double));
        // MORSE_Tile_to_Lapack( data->descC, C, N);
        // print_dmatrix("testC", N, N, C, N);
        //exit(0);

        //MORSE_MLE_dprint_Tile_Async(MORSE_descC, msequence, &mrequest[0]);
        // MORSE_Sequence_Wait(msequence);
        //exit(0);

        //exit(0);
        //Calculate Cholesky Factorization (C=LL-1)
        VERBOSE("Cholesky factorization of Sigma...");
        START_TIMING(time_facto);
        success = MORSE_dpotrf_Tile(MorseLower, MORSE_descC);
        STOP_TIMING(time_facto);
        SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
        flops = flops + FLOPS_DPOTRF(N);
        VERBOSE(" Done.\n");


        //double *C = (double *) malloc(N * N * sizeof(double));
        //MORSE_Tile_to_Lapack( MORSE_descC, C, N);
        //print_dmatrix("testC", 16, 16, C, 16);

        //*********************************************
        //you need to generate the full matrix
        /*MORSE_desc_t *MORSE_descC2       = NULL;
          MORSE_desc_t *MORSE_descC3       = NULL;
          MORSE_desc_t *MORSE_descC4       = NULL;
          MORSE_Sequence_Wait(msequence);
          EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC2,NULL , MorseRealDouble, 560, 560, 560 * 560, N, N, 0, 0, N, N, 1, 1);
          EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC4,NULL , MorseRealDouble, 560, 560, 560 * 560, N, N, 0, 0, N, N, 1, 1);
          MORSE_dlaset_Tile(MorseUpperLower, 0, 0, MORSE_descC4);
          MORSE_dlacpy_Tile(MorseLower ,MORSE_descC,MORSE_descC4);
          EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC3,NULL , MorseRealDouble, 560, 560, 560 * 560, N, N, 0, 0, N, N, 1, 1);
          MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC2, msequence, &mrequest[0], &data->l1, &data->l1, (double *)theta,  data->dm);
          MORSE_Sequence_Wait(msequence);
          MORSE_dgemm_Tile (MorseNoTrans, MorseTrans, 1, MORSE_descC4, MORSE_descC4, 0, MORSE_descC3);

          double error=0;
          double norm_c=0;
          MORSE_dgeadd_Tile_Async( MorseTrans, 1.0, MORSE_descC2
          ,-1.0, MORSE_descC3, msequence, &mrequest[0] );
          MORSE_Sequence_Wait(msequence);
          MORSE_dlange_Tile_Async( MorseFrobeniusNorm,
          MORSE_descC3, &error, msequence, &mrequest[0]
          );
          MORSE_Sequence_Wait(msequence);
          MORSE_dlange_Tile_Async( MorseFrobeniusNorm,
          MORSE_descC2, &norm_c, msequence, &mrequest[0]
          );


          MORSE_Sequence_Wait(msequence);
          printf("error: %e\n", (error/norm_c));
          exit(0);
          */
        //***************************************
        //MORSE_Tile_to_Lapack( MORSE_descC, C, N);
        //print_dmatrix("testC", 16, 16, C, 16);
        //exit(0);
        //MORSE_MLE_dprint_Tile_Async(MORSE_descC, msequence, &mrequest[0]);
        //MORSE_Sequence_Wait(msequence);
        //exit(0);

        //Calculate log(|C|) --> log(square(|L|))
        VERBOSE("Calculating the log determinant ...");
        START_TIMING(logdet_calculate);
        MORSE_MLE_dmdet_Tile_Async(MORSE_descC, msequence, &mrequest[0], MORSE_descdet);
        MORSE_Sequence_Wait(msequence);
        // printf("det: %f\n", data->det);
        logdet= 2*data->det;
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n");

        //        printf("logdet: %f\n",logdet);        


        // double *C = (double *) malloc(N * N * sizeof(double));
        // MORSE_Tile_to_Lapack( MORSE_descC, C, N);
        // print_dmatrix("testC", 70, 70, C, 70);
        // double *zz = MORSE_descZ->mat;
        // int k=0;
        //for(k=0;k<N;k++)
        //        printf("%f, %d\n", zz[k], k);
        //exit(0);

        //Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("Solving the linear system ...\n");
        START_TIMING(time_solve);
        MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ);
        STOP_TIMING(time_solve);
        flops = flops + FLOPS_DTRSM(MorseLeft,N, NRHS);
        VERBOSE(" Done.\n");    

        //double *zz = MORSE_descZ->mat;
        //	int k=0;
        //	for(k=0;k<N;k++)
        //	      printf("%f, %d\n", zz[k], k);
        //	  exit(0);


        //Calculate MLE likelihood
        VERBOSE("Calculating the MLE likelihood function ...");
        //dotp=0;
        //MORSE_MLE_core_ddotp_Async(MORSE_descZ,MORSE_descproduct,msequence, &mrequest[0]);
        //MORSE_Sequence_Wait(msequence);        

        MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct); 

        //***************************************

        if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0 )
        {

            loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ1, MORSE_descZ1, 0, MORSE_descproduct1);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ2, MORSE_descZ2, 0, MORSE_descproduct2);
            data->variance1 = (1.0/(N/2)) * data->dotp1;
            data->variance2 = (1.0/(N/2)) * data->dotp2;

            /*                loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
                              double *z = (double *) malloc(N * sizeof(double));
                              double *z1 = (double *) malloc((N/2) * sizeof(double));
                              double *z2 = (double *) malloc((N/2) * sizeof(double));
                              int dts=320;
                              int p_grid=1;
                              int q_grid=1;
                              MORSE_desc_t *MORSE_descZ1p;
                              MORSE_desc_t *MORSE_descZ2p;

                              MORSE_Tile_to_Lapack( MORSE_descZ, z, N);
                              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ1p, NULL, MorseRealDouble, dts, dts, dts * dts, N/2, 1, 0, 0, N/2, 1, p_grid, q_grid);
                              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ2p, NULL, MorseRealDouble, dts, dts, dts * dts, N/2, 1, 0, 0, N/2, 1, p_grid, q_grid);
                              int i=0;
                              int j=0;


                              for(i=0;i<N/2;i++)
                              {
                              z1[j]=z[i];
                              z2[j]=z[i+N/2];
            //      printf("%f- %f\n", z1[j], z2[j]);
            j++;
            }


            MORSE_Lapack_to_Tile( z1, N/2, MORSE_descZ1p);
            MORSE_Lapack_to_Tile( z2, N/2, MORSE_descZ2p);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ1p, MORSE_descZ1p, 0, MORSE_descproduct1);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ2p, MORSE_descZ2p, 0, MORSE_descproduct2);
            data->variance1 = (1.0/(N/2)) * data->dotp1;
            data->variance2 = (1.0/(N/2)) * data->dotp2;
            */
        }
        else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
        {

            loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);

            //to be optimized
            MORSE_stride_vec_Tile_Async(MORSE_descZ, MORSE_descZ1, MORSE_descZ2, msequence, &mrequest[0]);
            //MORSE_Sequence_Wait(msequence);
            //double *z = (double *) malloc(N * sizeof(double));
            //MORSE_Tile_to_Lapack( MORSE_descZ, z, N);
            //MORSE_Lapack_to_Tile( z, N/2, MORSE_descZ1);
            //MORSE_Lapack_to_Tile( &z[N/2], N/2, MORSE_descZ2);
            //*********************************
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ1, MORSE_descZ1, 0, MORSE_descproduct1);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ2, MORSE_descZ2, 0, MORSE_descproduct2);
            data->variance1 = (1.0/(N/2)) * data->dotp1;
            data->variance2 = (1.0/(N/2)) * data->dotp2;
            //free(z);
            /*                loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
                              double *z = (double *) malloc(N * sizeof(double));
                              double *z1 = (double *) malloc((N/2) * sizeof(double));
                              double *z2 = (double *) malloc((N/2) * sizeof(double));
                              int dts=320;
                              int p_grid=1;
                              int q_grid=1;
                              MORSE_desc_t *MORSE_descZ1p;
                              MORSE_desc_t *MORSE_descZ2p;

                              MORSE_Tile_to_Lapack( MORSE_descZ, z, N);
                              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ1p, NULL, MorseRealDouble, dts, dts, dts * dts, N/2, 1, 0, 0, N/2, 1, p_grid, q_grid);
                              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ2p, NULL, MorseRealDouble, dts, dts, dts * dts, N/2, 1, 0, 0, N/2, 1, p_grid, q_grid);
                              int i=0;
                              int j=0;


                              for(i=0;i<N;i+=2)
                              {
                              z1[j]=z[i];
                              z2[j]=z[i+1];
            //	printf("%f- %f\n", z1[j], z2[j]);
            j++;
            }


            MORSE_Lapack_to_Tile( z1, N/2, MORSE_descZ1p);
            MORSE_Lapack_to_Tile( z2, N/2, MORSE_descZ2p);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ1p, MORSE_descZ1p, 0, MORSE_descproduct1);
            MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ2p, MORSE_descZ2p, 0, MORSE_descproduct2);
            data->variance1 = (1.0/(N/2)) * data->dotp1;
            data->variance2 = (1.0/(N/2)) * data->dotp2;
            */

        }
        else
        {

            loglik = -0.5 * data->dotp -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
            data->variance= theta[0];
        }
        VERBOSE(" Done.\n");

    }

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    //MPI_Bcast(theta, num_params, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(MORSE_My_Mpi_Rank() == 0)
    {
#endif

        //		if(strcmp(data->checkpoint_file,"") != 0)    
        //			checkpointing(data->checkpoint_file, data->iter_count, theta, loglik, num_params);


        //Print Iteration Summary
        //fprintf(stderr,"***************************************************\n");
        //	fprintf(stderr,"\n------ddotproduct: %.17g ", data->dotp);
        //	fprintf(stderr,"\n------logdet: %.17g ", logdet);
        //fprintf(stderr,"------det: %.*e ", det);
        //fprintf(stderr,"\n------expr2: %.8f \n",((double) (N / 2) * log(2 * PI)));
        //fprintf(stderr," ---- Theta1: %.8f ----  Theta2: %.8f ---- Theta3: %.8f ----LogLi: %.8f\n", theta[0], theta[1], theta[2],loglik);
        //reformat

        printf(" %3d- Model Parameters (",  data->iter_count+1);

        if(data->log == 1)
            fprintf(data->pFileLog, " %3d- Model Parameters (",  data->iter_count+1);

        if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        {	
            printf("%.8f, %.8f,", data->variance1, data->variance2);
            if(data->log == 1)
                fprintf(data->pFileLog,"%.8f, %.8f,", data->variance1, data->variance2);
            i = 2;
            results.estimated_theta[0] = data->variance1;
            results.estimated_theta[1] = data->variance2;
        }
        else
            i=0;
        for(;i<num_params; i++)
        {
            printf("%.8f", theta[i]);
            if (i <num_params-1)
                printf(",");

            results.estimated_theta[i] = theta[i];
            if(data->log == 1)
                fprintf(data->pFileLog,"%.8f, ", theta[i]);
        }

        printf(")----> LogLi: %.18f\n", loglik);
        if(data->log == 1)
            fprintf(data->pFileLog, ")----> LogLi: %.18f\n", loglik);


        printf(" ---- Facto Time: %6.2f\n", time_facto);
        printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
        printf(" ---- dtrsm Time: %6.2f\n", time_solve);
        printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
        //fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", zcpy_time);
        printf(" ---- Total Time: %6.2f\n", /*matrix_gen_time+*/ time_facto + logdet_calculate + time_solve);
        //fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );
        printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
        //fprintf(stderr," ---- Peak Performance: %6.2f Gflops/s\n",  (ncores*p_grid*q_grid*16*2.3) );
        //fprintf(stderr,"***************************************************\n");

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter += /*matrix_gen_time +*/ time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    data->final_loglik = loglik;

    //output
    results.final_loglik = loglik;


    return loglik;
}

double MORSE_dmle_Tile_Async(unsigned n, const double * theta, double * grad, void * MORSE_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- MORSE-Async
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library.
     * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0, dzcpy_time=0.0, flops = 0.0;
    int N, NRHS, success;

    MLE_data* data	= ((MLE_data*)MORSE_data);
    data->det	= 0;
    data->dotp	= 0;

    MORSE_desc_t *MORSE_descC	= (MORSE_desc_t *) data->descC;
    MORSE_desc_t *MORSE_descZ	= (MORSE_desc_t *) data->descZ;
    MORSE_desc_t *MORSE_descZcpy	= (MORSE_desc_t *) data->descZcpy; 
    MORSE_desc_t *MORSE_descdet	= (MORSE_desc_t *) data->descdet;
    MORSE_desc_t *MORSE_descproduct	= (MORSE_desc_t *) data->descproduct;
    MORSE_sequence_t *msequence	= (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest	= (MORSE_request_t *) data->request;

    N	= MORSE_descC->m;
    NRHS	= MORSE_descZ->n;
    START_TIMING(dzcpy_time);
    if(data->iter_count == 0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descZ, MORSE_descZcpy, msequence, mrequest); 
    else
    {	
        VERBOSE("re-store the original Z vector...");
        MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descZcpy, MORSE_descZ, msequence, mrequest);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(dzcpy_time);	


    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix...");
    START_TIMING(matrix_gen_time);	
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest, &data->l1, &data->l1,(double*) theta,  data->dm);    
    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, &data->l1, &data->l1, &data->lm, (double *)theta, data->dm, data->kernel_fun,  msequence, &mrequest[0]);
    MORSE_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma...");
    START_TIMING(time_facto);
    success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_DPOTRF(N);
    VERBOSE(" Done.\n");

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant ...");
    START_TIMING(logdet_calculate);
    MORSE_MLE_dmdet_Tile_Async(MORSE_descC, msequence, &mrequest[0],MORSE_descdet);
    MORSE_Sequence_Wait(msequence);
    logdet= 2*data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system ...\n");
    START_TIMING(time_solve);
    MORSE_dtrsm_Tile_Async(MorseLeft,MorseLower,MorseNoTrans,MorseNonUnit,1,MORSE_descC,MORSE_descZ, msequence, mrequest);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_DTRSM(MorseLeft,N, NRHS);
    VERBOSE(" Done.\n");    

    //Claculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...");
    MORSE_dgemm_Tile_Async (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct, msequence, mrequest); 
    loglik = -0.5 * data->dotp -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
    VERBOSE(" Done.\n");

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(MORSE_My_Mpi_Rank() == 0)
    {
#endif
        //Print Iteration Summary
        //fprintf(stderr,"***************************************************\n");
        fprintf(stderr,"------ddotproduct: %.8f ", data->dotp);
        fprintf(stderr,"------logdet: %.8f ", logdet);
        //fprintf(stderr,"------det: %.*e ", det);
        fprintf(stderr,"------expr2: %.8f ",((double) (N / 2) * log(2 * PI)));
        //fprintf(stderr," ---- Theta1: %.8f ----  Theta2: %.8f ---- Theta3: %.8f ----LogLi: %.8f\n", theta[0], theta[1], theta[2],loglik);
        //reformat
        fprintf(stderr," %3d- Model Parameters (variance, range, smoothness): (%.8f, %.8f, %.8f) ----> LogLi: %.8f\n", data->iter_count+1,  theta[0], theta[1], theta[2],loglik);

        if(data->log == 1)
            fprintf(data->pFileLog, " %3d- Model Parameters (variance, range, smoothness): (%.8f, %.8f, %.8f) ----> LogLi: %.8f\n", data->iter_count+1,  theta[0], theta[1], theta[2],loglik);

        fprintf(stderr," ---- Facto Time: %6.2f\n", time_facto);
        fprintf(stderr," ---- logdet Time: %6.2f\n", logdet_calculate);
        fprintf(stderr," ---- dtrsm Time: %6.2f\n", time_solve);
        fprintf(stderr," ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
        //fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", dzcpy_time);
        fprintf(stderr," ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
        //fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );    
        //fprintf(stderr," ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
        //fprintf(stderr," ---- Peak Performance: %6.2f Gflops/s\n",  (ncores*p_grid*q_grid*16*2.3) );
        //fprintf(stderr,"***************************************************\n");
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter+=matrix_gen_time+time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter+=flops / 1e9 / (time_facto +time_solve);
    data->final_loglik=loglik;

    return loglik;
}


void MORSE_dmle_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid, int mse_flag)
    //! Allocate prediction operation descriptors.
    /*!  
     * Returns MLE_data data with initial values and new descriptors locations.
     * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
     * @param[in] nZmiss: number of missing values (unknown observations).
     * @param[in] nZobs: number of observed values (known observations).
     * @param[in] dts: tile size (MB).
     * @param[in] p_grid: p_grid in the case of distributed system.
     * @param[in] q_grid: q_grid in the case of distributed system.
     * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
     * */
{

    MORSE_desc_t *MORSE_descZmiss   = NULL;
    MORSE_desc_t *MORSE_descC12     = NULL;
    MORSE_desc_t *MORSE_descC22     = NULL;
    MORSE_desc_t *MORSE_descmse     = NULL;
    MORSE_desc_t *MORSE_descmse1     = NULL;
    MORSE_desc_t *MORSE_descmse2     = NULL;	
    MORSE_desc_t *MORSE_descZactual = NULL;
    MORSE_desc_t *MORSE_descZobs    = NULL;
    MLE_data     *data              = (MLE_data*) MORSE_data;

    if(nZmiss <= 0)
    {
        fprintf(stderr," Number of missing values should be positive value\n");
        return;
    }

    //bi-variate case
    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_flexible_profile")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_flexible2")   == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_flexible2_profile")   == 0 )
    {

        nZobs*=2;
        nZmiss*=2;

    }


    //Descriptors Creation
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
    if( mse_flag == 1)
    {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse1, &data->mserror1, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse2, &data->mserror2, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

    //Initiate data descriptors
    data->descZmiss         = MORSE_descZmiss;
    data->descC12           = MORSE_descC12;
    data->descC22           = MORSE_descC22;
    data->descmse           = MORSE_descmse;
    data->descmse1           = MORSE_descmse1;
    data->descmse2           = MORSE_descmse2;
    data->descZactual       = MORSE_descZactual;
    data->descZobs          = MORSE_descZobs;

    //printf("%d- %d - %d\n", MORSE_descC12->m, MORSE_descC12->n, MORSE_descZobs->m);
    //exit(0);

}


double MORSE_dmle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n)
    //! //Predict missing values base on a set of given values and covariance matrix
    /*!  -- MORSE-sync
     * Returns the prediction Mean Square Error (MSE) as double
     * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
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
    //double *z = NULL, *streamdata = NULL;
    double time_solve = 0.0;
    double mat_gen_time = 0.0;
    double time_gemm = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;
    int num_params = 0;

    MORSE_desc_t *MORSE_descZmiss   = NULL;
    MORSE_desc_t *MORSE_descC12     = NULL;
    MORSE_desc_t *MORSE_descC22     = NULL;
    MORSE_desc_t *MORSE_descmse     = NULL;
    MORSE_desc_t *MORSE_descmse1     = NULL;
    MORSE_desc_t *MORSE_descmse2     = NULL;	
    MORSE_desc_t *MORSE_descZactual = NULL;
    MORSE_desc_t *MORSE_descZobs    = NULL;
    MLE_data     *data              = (MLE_data*) MORSE_data;
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t *mrequest       = (MORSE_request_t  *) data->request;
    data->mserror                   = 0;

    if(nZmiss <= 0)
    {
        fprintf(stderr," Number of missing values should be positive value\n");
        return -1;
    }

    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
    {
        data->kernel_fun= "bivariate_matern_parsimonious";
    }
    //Descriptors Creation
    //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
    //if( Zactual != NULL)
    //{
    //        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
    //        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    //}
    //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
    //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

    //Initiate data descriptors
    MORSE_descZmiss		= data->descZmiss;
    MORSE_descC12		= data->descC12;
    MORSE_descC22		= data->descC22;
    MORSE_descmse		= data->descmse;
    MORSE_descmse1           = data->descmse1;
    MORSE_descmse2           = data->descmse2;
    MORSE_descZactual	= data->descZactual;
    MORSE_descZobs		= data->descZobs;

    //Copy data to vectors 
    VERBOSE("Copy measurments vector to descZobs descriptor...");
    //MORSE_MLE_dzcpy_Tile_Async(MORSE_descZobs, Zobs, msequence, mrequest);
    MORSE_Lapack_to_Tile( Zobs, nZobs, MORSE_descZobs);
    VERBOSE(" Done.\n");

    if( Zactual != NULL)	
    {
        VERBOSE("Copy actual measurments vector to descZactual descriptor...");
        //MORSE_MLE_dzcpy_Tile_Async(MORSE_descZactual, Zactual, msequence, mrequest);
        MORSE_Lapack_to_Tile( Zactual, nZmiss, MORSE_descZactual);
        VERBOSE(" Done.\n");
    }

    MORSE_Sequence_Wait(msequence);


    //        int i=0;
    //      for (i=0;i<nZmiss;i++)
    //    printf("%f, %f, %f\n", data->lmiss.x[i], data->lmiss.y[i], Zactual[i]);

    //printf("\n\n");

    //	for (i=0;i<100;i++)
    //	printf("%f, %f, %f\n", data->lobs.x[i], data->lobs.y[i], Zobs[i]);

    //#if defined(CHAMELEON_USE_MPI)
    //	MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    //#endif
    //theta[0]=data->variance;
    if(strcmp(data->kernel_fun, "univariate_matern_stationary")   == 0 || strcmp(data->kernel_fun, "univariate_pow_exp_stationary")   == 0 )
        num_params = 3;
    else if(strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary")   == 0)
        num_params = 4;
    else if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
        num_params = 9;
    else if(strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0)
        num_params = 11;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        num_params = 6;
    else if(strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary")   == 0)
        num_params = 7;
    else
    {
        fprintf(stderr,"Choosen kernel is not exist(1)!\n");
        fprintf(stderr, "Called function is: %s\n",__func__);
        exit(0);
    }

    printf("estimated parameters:");
    int i = 0;
    for(i=0; i<num_params; i++)
    {
        printf("%.8f,", theta[i]);

    }
    printf(")\n");
    START_TIMING(mat_gen_time);




    //double *Zobs_arr = (double *) malloc(nZobs * 1 * sizeof(double));
    //MORSE_Tile_to_Lapack( MORSE_descZobs, Zobs_arr, nZobs);
    //print_dmatrix("Zobs(16)", 16, 1, Zobs_arr, 1);

    //double *Zactual_arr = (double *) malloc(nZmiss * 1 * sizeof(double));
    //MORSE_Tile_to_Lapack( MORSE_descZactual, Zactual_arr, nZmiss);
    //print_dmatrix("Zactual(10)", 10, 1, Zactual_arr, 1);



    //Generate C22 covariance matrix
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  &data->lobs, &data->lobs, theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZobs);
    VERBOSE(" Done.\n");


    //printf("%s\n",  data->kernel_fun);
    //double *C22_arr = (double *) malloc(nZobs * nZobs * sizeof(double));
    //MORSE_Tile_to_Lapack( MORSE_descC22, C22_arr, nZobs);
    //print_dmatrix("C22(10,10)", 10, 10, C22_arr, 10);

    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  &data->lmiss, &data->lobs, theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, &data->lmiss, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZmiss);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);


    //printf("%s\n",  data->kernel_fun);
    //double *C12_arr = (double *) malloc(nZobs * nZmiss * sizeof(double));
    //MORSE_Tile_to_Lapack( MORSE_descC12, C12_arr, nZmiss);
    //print_dmatrix("C12(10,10)", 10, 10, C12_arr, 10);

    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
    MORSE_dposv_Tile(MorseLower, MORSE_descC22, MORSE_descZobs);
    flops = flops + FLOPS_DPOTRF(nZobs);
    flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descC22->m, MORSE_descZobs->n);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);


    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
    MORSE_dgemm_Tile (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss);
    flops = flops + FLOPS_DGEMM(MORSE_descC12->m, MORSE_descZobs->n, MORSE_descZmiss->n);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //printf("%s\n",  data->kernel_fun);
    //double *Zmiss_arr = (double *) malloc(nZmiss * 1 * sizeof(double));
    //MORSE_Tile_to_Lapack( MORSE_descZmiss, Zmiss_arr, nZobs);
    //print_dmatrix("Zmiss(10)", 10, 1, Zmiss_arr, 1);
    //	exit(0);

    //return back descZmiss to zmiss vector
    MORSE_Tile_to_Lapack( MORSE_descZmiss, Zmiss, nZmiss);

    //Estimate Mean Square Error
    if( Zactual != NULL)
    {
	    START_TIMING(time_mse);
	    VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");

	    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0 ) 

		    MORSE_MLE_dmse_bivariate_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse1, MORSE_descmse2, MORSE_descmse,  msequence, mrequest);
	    else

		    MORSE_MLE_dmse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
	    MORSE_Sequence_Wait(msequence);
	    VERBOSE(" Done.\n");	
	    STOP_TIMING(time_mse);
	    data->mserror  /= nZmiss;
	    data->mserror1 /= (nZmiss/2);
	    data->mserror2 /= (nZmiss/2);
    }
    else
	    data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
    if(MORSE_My_Mpi_Rank() == 0)
    {
#endif

	    double *z;
	    //		z = (double *) malloc(nZmiss * sizeof(double));
	    //		MORSE_Tile_to_Lapack( MORSE_descZmiss, z, nZmiss);
	    //		write_pred_vector(z, data, nZmiss, data->mserror);
	    //		free(z);
	    if(data->log == 1)
		    fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, Flops: %.8f, Mean Square Error (MSE): %.8f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

	    write_prediction_result("predict_result.dat", n, data->hicma_acc, data->mserror1, data->mserror2, data->mserror, (mat_gen_time+time_solve+ time_gemm), (flops / 1e9 / (time_solve)));


	    //output
	    results.mse_pred1 = data->mserror1;
	    results.mse_pred2 = data->mserror2;
	    results.mse_pred  = data->mserror;
	    results.total_pred_time= /*mat_gen_time+*/time_solve+ time_gemm;
	    results.total_pred_flops= flops / 1e9 / (time_solve + time_gemm);

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;

}



double MORSE_dmle_Predict_Tile_Async(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-Async
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
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
	MLE_data * data			=  (MLE_data*)MORSE_data;
	MORSE_desc_t *MORSE_descZ    	= (MORSE_desc_t *)(data->descZcpy);
	MORSE_desc_t *MORSE_descZobs	= (MORSE_desc_t *)(data->descZobs);	
	MORSE_desc_t *MORSE_descZactual	= (MORSE_desc_t *)(data->descZactual);
	MORSE_desc_t *MORSE_descZmiss 	= (MORSE_desc_t *)(data->descZmiss);	
	MORSE_desc_t *MORSE_descC12 	= (MORSE_desc_t *)(data->descC12);
	MORSE_desc_t *MORSE_descC22 	= (MORSE_desc_t *)(data->descC22);
	MORSE_desc_t *MORSE_descmse 	= (MORSE_desc_t *)(data->descmse);
	MORSE_sequence_t *msequence 	= (MORSE_sequence_t *)(data->sequence);
	MORSE_request_t *mrequest 	= (MORSE_request_t *)data->request;

	if(strcmp(data->actualZFPath,"")==0)
	{
		double *z = NULL;
#if defined(CHAMELEON_USE_MPI)
		z = (double *) malloc(n * sizeof(double));
		MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
#else
		z = MORSE_descZ->mat;
#endif

		//random  shuffle
		//	shuffle(z, &data->l1, n);

#if defined(CHAMELEON_USE_MPI)
		MORSE_Lapack_to_Tile( z, n, MORSE_descZ);
#endif

		l1 = &data->l1;
		temp_loc.x=&l1->x[nZmiss];
		temp_loc.y=&l1->y[nZmiss];
		l2 = &temp_loc;
	}
	else
	{
		double *streamdata = NULL;
		l1 = &data->l1;
		temp_loc.x=&l1->x[nZmiss];
		temp_loc.y=&l1->y[nZmiss];
		l2 = &temp_loc;

		//l1 = (location *) malloc(sizeof(location));

		//l1->x=(double *) malloc(nZmiss * sizeof(double));
		//l1->y=(double *) malloc(nZmiss * sizeof(double));

		VERBOSE("Reading ActualZ locations for prediction from disk .....");
		l1 = readLocsFile(data->actualZLocFPath, nZmiss);	
		VERBOSE(" Done.\n");

		//streamdata=(double *) malloc(nZmiss * sizeof(double));
		VERBOSE("Reading ActualZ for prediction from disk .....");
		streamdata = readObsFile(data->actualZFPath, nZmiss);
		MORSE_MLE_dzcpy_Tile_Async(MORSE_descZactual, streamdata, msequence, mrequest);
		MORSE_Sequence_Wait(data->sequence);
		VERBOSE(" Done.\n");
	}


	//MORSE_dposv_Tile_Async(MorseLower, MORSE_descC22, MORSE_descZobs, data->sequence, &data->request[0]);
	START_TIMING(mat_gen_time);

	//Generate C22 covariance matrix
	VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
	//MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  l2, l2, theta, data->dm);
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC22, l2, l2, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
	//flops = flops + FLOPS_DPOTRF(nZobs);
	VERBOSE(" Done.\n");

	//Generate C12 covariance matrix
	VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
	//MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  l1, l2, theta, data->dm);
	MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, l1, l2, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
	//flops = flops + FLOPS_DPOTRF(nZmiss);
	VERBOSE(" Done.\n");
	STOP_TIMING(mat_gen_time);

	START_TIMING(time_solve);
	//Start prediction
	VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
	MORSE_dposv_Tile_Async(MorseLower, MORSE_descC22, MORSE_descZobs, msequence, mrequest);
	flops = flops + FLOPS_DPOTRF(nZobs);
	flops = flops + FLOPS_DTRSM(MorseLeft, nZobs, nZobs);
	VERBOSE(" Done.\n");

	VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
	MORSE_dgemm_Tile_Async (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss, msequence, mrequest);
	flops = flops + FLOPS_DGEMM(MORSE_descC12->m, MORSE_descZobs->n, MORSE_descZmiss->n);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_solve);


	//Estimate Mean Square Error
	START_TIMING(time_mse);
	VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
	MORSE_MLE_dmse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
	VERBOSE(" Done.\n");
	STOP_TIMING(time_mse);


	//if you do not have actual value to compare with
	if(data->descZactual==NULL)
		return -1;

	data->mserror /= nZmiss;

#if defined(CHAMELEON_USE_MPI)
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		if(data->log == 1)
			fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %.8f, Flops: %.8f, Mean Square Error (MSE): %.8f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

		// write_prediction_result("predict_result.dat", n, nZmiss, data->mserror, (mat_gen_time+time_solve+ time_mse), (flops / 1e9 / (time_solve )));

#if defined(CHAMELEON_USE_MPI)
	}
#endif

	return data->mserror;
}



void MORSE_dmle_mloe_mmom_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid)
	//! Allocate prediction operation descriptors.
	/*!
	 * Returns MLE_data data with initial values and new descriptors locations.
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
	 * @param[in] nZmiss: number of missing values (unknown observations).
	 * @param[in] nZobs: number of observed values (known observations).
	 * @param[in] dts: tile size (MB).
	 * @param[in] p_grid: p_grid in the case of distributed system.
	 * @param[in] q_grid: q_grid in the case of distributed system.
	 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
	 * */
{

	MORSE_desc_t *MORSE_desck_t        	= NULL;
	MORSE_desc_t *MORSE_desck_a         	= NULL;
	MORSE_desc_t *MORSE_desck_atmp      	= NULL;
	MORSE_desc_t *MORSE_desck_ttmp      	= NULL;
	//MORSE_desc_t *MORSE_desck_atmp3     = NULL;

	MORSE_desc_t *MORSE_descK_t         	= NULL;
	MORSE_desc_t *MORSE_descK_ttmp      	= NULL;
	MORSE_desc_t *MORSE_descK_a         	= NULL;

	MORSE_desc_t *MORSE_descexpr1       	= NULL;
	MORSE_desc_t *MORSE_descexpr2       	= NULL;
	MORSE_desc_t *MORSE_descexpr3       	= NULL;
	MORSE_desc_t *MORSE_descexpr4       	= NULL;
	MORSE_desc_t *MORSE_descalpha       	= NULL;
	MORSE_desc_t *MORSE_desctruthalpha      = NULL;
	MORSE_desc_t *MORSE_descestimatedalpha	= NULL;
	MORSE_desc_t *MORSE_desc_mloe_mmom      = NULL;
	MLE_data     *data              = (MLE_data*) MORSE_data;
	int p=0;
	if(nZmiss <= 0)
	{
		fprintf(stderr," Number of missing values should be positive value\n");
		return;
	}
	if(strcmp(data->kernel_fun, "univariate_matern_stationary")   == 0)
		p	= 1;
	else if(strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary")   == 0)
		p       = 1;
	else if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
		p       = 1;
	else if(strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0)
		p       = 2;
	else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
		p       = 2;
	else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
		p       = 2;
	else if(strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary")   == 0)
		p	= 1;
	else
	{
		fprintf(stderr,"Choosen kernel is not exist(24)!\n");
		fprintf(stderr, "Called function is: %s\n",__func__);
		exit(0);
	}



	//EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desc_mloe_mmom, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_t,    NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p,  0, 0, p*nZobs, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_a,    NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p,  0, 0, p*nZobs, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_atmp, NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p,  0, 0, p*nZobs, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desck_ttmp, NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p,  0, 0, p*nZobs, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr1, NULL, MorseRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr2, NULL, MorseRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr3, NULL, MorseRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descexpr4, NULL, MorseRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_desctruthalpha, NULL, MorseRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descestimatedalpha, NULL, MorseRealDouble, dts, dts, dts * dts, p, p, 0, 0, p, p, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descK_t, NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p*nZobs, 0, 0, p*nZobs, p*nZobs, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descK_a, NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p*nZobs, 0, 0, p*nZobs, p*nZobs, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descK_ttmp, NULL, MorseRealDouble, dts, dts, dts * dts, p*nZobs, p*nZobs, 0, 0, p*nZobs, p*nZobs, p_grid, q_grid);

	//Initiae data descriptors
	data->desck_t			= MORSE_desck_t;
	data->desck_a			= MORSE_desck_a;
	data->descK_ttmp		= MORSE_descK_ttmp;
	data->desck_atmp       		= MORSE_desck_atmp;
	data->desck_ttmp        	= MORSE_desck_ttmp;
	data->descK_t			= MORSE_descK_t;
	data->descK_a			= MORSE_descK_a;
	data->descexpr1			= MORSE_descexpr1;
	data->descexpr2			= MORSE_descexpr2;
	data->descexpr3			= MORSE_descexpr3;
	data->descexpr4			= MORSE_descexpr4;
	data->descestimatedalpha        = MORSE_descestimatedalpha;
	data->desctruthalpha        	= MORSE_desctruthalpha;
	//data->desc_mloe_mmom    	= MORSE_desc_mloe_mmom;
}


void MORSE_dmle_mloe_mmom_Tile(MLE_data *MORSE_data, double * truth_theta, double* estimated_theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-sync
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
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
	//truth_theta[0]=1; truth_theta[1]=0.1; truth_theta[2]=0.5;
	//estimated_theta[0]=1.01; estimated_theta[1]=0.09; estimated_theta[2]=0.49;
	printf("%f, %f, %f, %f\n", truth_theta[0], truth_theta[1], truth_theta[2], truth_theta[3]);
	printf("%f, %f, %f,%f \n", estimated_theta[0], estimated_theta[1], estimated_theta[2], estimated_theta[3]);
	//exit(0);     
	double loe_sum	= 0.0;
	double mom_sum	= 0.0;
	int i		= 0;
	int p		= 0;
	int v 		= 0;
	int j 		= 0;	
	double all_time = 0.0;
	double cholesky1 = 0.0;
	double cholesky2 = 0.0;
	double matrix_gen1 = 0.0;
	double matrix_gen2 = 0.0;
	double matrix_gen3 = 0.0;
	double matrix_gen4 = 0.0;
	double copy1 = 0.0;
	double copy2 = 0.0;
	double copy3 = 0.0;

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
	double *loe=(double *) malloc(nZmiss * sizeof(double));
	double *mom=(double *) malloc(nZmiss * sizeof(double));

	MLE_data     *data                	  = (MLE_data*) MORSE_data;
	MORSE_desc_t *MORSE_desck_t        	  = data->desck_t;
	MORSE_desc_t *MORSE_desck_a         	  = data->desck_a;
	MORSE_desc_t *MORSE_descK_t         	  = data->descK_t;
	MORSE_desc_t *MORSE_descK_a               = data->descK_a;
	MORSE_desc_t *MORSE_descK_ttmp            = data->descK_ttmp;
	MORSE_desc_t *MORSE_desck_atmp            = data->desck_atmp;
	MORSE_desc_t *MORSE_desck_ttmp            = data->desck_ttmp;
	MORSE_desc_t *MORSE_descexpr1             = data->descexpr1;
	MORSE_desc_t *MORSE_descexpr2             = data->descexpr2;
	MORSE_desc_t *MORSE_descexpr3             = data->descexpr3;
	MORSE_desc_t *MORSE_descexpr4             = data->descexpr4;
	MORSE_desc_t *MORSE_descestimatedalpha    = data->descestimatedalpha;
	MORSE_desc_t *MORSE_desctruthalpha        = data->desctruthalpha;
	MORSE_sequence_t *msequence  	          = (MORSE_sequence_t *)(data->sequence);
	MORSE_request_t *mrequest                 = (MORSE_request_t *)data->request;
	location lmiss; 
	lmiss.x            = (double *) malloc(sizeof(double));
	lmiss.y            = (double *) malloc(sizeof(double));

	double* univariate_theta;
	double* univariate2_theta;
	double* univariate3_theta;
	double nu12;
	double rho;
	double sigma_square12;

	double flops = 0.0;
	START_TIMING(all_time);

	int m = MORSE_descestimatedalpha->m;
	double *truthalpha    = (double *) malloc( m * m * sizeof(double));
	double *estimatedalpha = (double *) malloc( m * m * sizeof(double));
	double *temp1 = (double *) malloc( m * m* sizeof(double));
	double *temp2 = (double *) malloc( m * m*sizeof(double));
	double *temp3 = (double *) malloc( m * m*sizeof(double));
	if(m ==1)
	{
		truthalpha[0]   = truth_theta[0];
		estimatedalpha[0] = estimated_theta[0];	                
	}

	if(m ==2)
	{
		double truth_nu12 = 0.5 * (truth_theta[3] + truth_theta[4]);
		double truth_rho = truth_theta[5] * sqrt( (tgamma(truth_theta[3] + 1)*tgamma(truth_theta[4] + 1)) /
				(tgamma(truth_theta[3]) * tgamma(truth_theta[4])) ) *
			tgamma(truth_nu12) / tgamma(truth_nu12 + 1);

		double estimated_nu12 = 0.5 * (estimated_theta[3] + estimated_theta[4]);
		double estimated_rho = estimated_theta[5] * sqrt( (tgamma(estimated_theta[3] + 1)*tgamma(estimated_theta[4] + 1)) /
				(tgamma(estimated_theta[3]) * tgamma(estimated_theta[4])) ) *
			tgamma(estimated_nu12) / tgamma(estimated_nu12 + 1);

		truthalpha[0]    = truth_theta[0];
		estimatedalpha[0] = estimated_theta[0];

		truthalpha[1]    = truthalpha[3] = truth_rho
			* sqrt(truth_theta[0] * truth_theta[1]);

		estimatedalpha[1] = estimatedalpha[3] = estimated_rho
			* sqrt(estimated_theta[0] * estimated_theta[1]);
		truthalpha[2]    = truth_theta[1];
		estimatedalpha[2] = estimated_theta[1];

	}
	MORSE_Lapack_to_Tile( truthalpha, m, MORSE_desctruthalpha);	
	MORSE_Lapack_to_Tile( estimatedalpha, m, MORSE_descestimatedalpha);

	char * name= malloc(strlen(data->kernel_fun) + 1);
	strcpy(name, data->kernel_fun);

	if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
	{
		printf("TODO:  running with Z=1 no more\n");
		data->kernel_fun= "bivariate_matern_parsimonious";
	}

	START_TIMING(matrix_gen1);
	VERBOSE("Create MORSE_descK_a Covariance Matrix (MLOE-MMOM).....");
	if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
	{

		univariate_theta =(double *) malloc(3 * sizeof(double));
		univariate2_theta =(double *) malloc(3 * sizeof(double));
		univariate3_theta =(double *) malloc(3 * sizeof(double));
		univariate_theta[0]=estimated_theta[0];
		univariate_theta[1]=estimated_theta[2];
		univariate_theta[2]=estimated_theta[3];

		MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, morse_desc_submatrix(MORSE_descK_a, 0,   0, MORSE_descK_a->m/2, MORSE_descK_a->n/2), &data->lobs, &data->lobs, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);


		nu12 = 0.5 * (estimated_theta[3] + estimated_theta[4]);
		rho = estimated_theta[5] * sqrt( (tgamma(estimated_theta[3] + 1)*tgamma(estimated_theta[4] + 1)) /
				(tgamma(estimated_theta[3]) * tgamma(estimated_theta[4])) ) *
			tgamma(nu12) / tgamma(nu12 + 1);
		sigma_square12 = rho * sqrt(estimated_theta[0]*estimated_theta[1]) ;

		univariate2_theta[0]=sigma_square12;
		univariate2_theta[1]=estimated_theta[2];
		univariate2_theta[2]=nu12;

		MORSE_MLE_dcmg_Tile_Async(MorseUpperLower,  morse_desc_submatrix(MORSE_descK_a, MORSE_descK_a->m/2,   0, MORSE_descK_a->m/2, MORSE_descK_a->n/2), &data->lobs, &data->lobs, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);
		//      MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descsubC21, &data->l1, &data->l1, &data->lm, univariate2_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

		univariate3_theta[0]=estimated_theta[1];
		univariate3_theta[1]=estimated_theta[2];
		univariate3_theta[2]=estimated_theta[4];



		MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, morse_desc_submatrix(MORSE_descK_a, MORSE_descK_a->m/2,   MORSE_descK_a->n/2, MORSE_descK_a->m/2, MORSE_descK_a->n/2), &data->lobs, &data->lobs, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);
	}
	else
	{		//(1)Generate the co-variance matrix descK_a
		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descK_a,  &data->lobs, &data->lobs, &data->lm, estimated_theta, data->dm, data->kernel_fun,  msequence, mrequest);
	}
	VERBOSE(" Done.\n");
	MORSE_Sequence_Wait(msequence);
	STOP_TIMING(matrix_gen1);	

	START_TIMING(matrix_gen2);
	//(2)Generate the co-variance matrix descK_t
	VERBOSE("Create MORSE_descK_t Covariance Matrix (MLOE-MMOM).....");
	MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descK_t,  &data->lobs, &data->lobs, &data->lm, truth_theta, data->dm, data->kernel_fun,  msequence, mrequest);
	MORSE_Sequence_Wait(msequence);
	VERBOSE(" Done.\n");
	STOP_TIMING(matrix_gen2);


	START_TIMING(cholesky1);
	//(3)Cholesky factorization for the Co-variance matrix MORSE_descK_a
	VERBOSE("Cholesky factorization of MORSE_descK_a (MLOE-MMOM) .....");
	int success = MORSE_dpotrf_Tile(MorseLower, MORSE_descK_a);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");
	STOP_TIMING(cholesky1);
	flops = flops + FLOPS_DPOTRF(MORSE_descK_a->m);


	START_TIMING(copy1);
	//(4)Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM)
	VERBOSE("Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM).....");
	MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descK_t, MORSE_descK_ttmp);
	VERBOSE(" Done.\n");
	STOP_TIMING(copy1);

	START_TIMING(cholesky2);
	//(5)Cholesky factorization for the Co-variance matrix MORSE_descK_t
	VERBOSE("Cholesky factorization of MORSE_descK_t (MLOE-MMOM) .....");
	success = MORSE_dpotrf_Tile(MorseLower, MORSE_descK_t);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");
	STOP_TIMING(cholesky2);
	flops = flops + FLOPS_DPOTRF(MORSE_descK_t->m);



	MORSE_desc_t *B2;
	MORSE_desc_t *B3;
	MORSE_desc_t *B4;
	MORSE_Desc_Create_User( &B2, temp1, MorseComplexDouble, MORSE_descexpr2->mb, MORSE_descexpr2->nb, MORSE_descexpr2->bsiz,
			m, MORSE_descexpr2->n, 0, 0, MORSE_descexpr2->m, MORSE_descexpr2->n, 1, 1,
			morse_getaddr_cm, morse_getblkldd_cm, NULL );
	MORSE_Sequence_Wait(msequence);	
	MORSE_Desc_Create_User( &B3, temp2, MorseComplexDouble, MORSE_descexpr3->mb, MORSE_descexpr3->nb, MORSE_descexpr3->bsiz,
			m, MORSE_descexpr3->n, 0, 0, MORSE_descexpr3->m, MORSE_descexpr3->n, 1, 1,
			morse_getaddr_cm, morse_getblkldd_cm, NULL );
	MORSE_Sequence_Wait(msequence);
	MORSE_Desc_Create_User( &B4, temp3, MorseComplexDouble, MORSE_descexpr4->mb, MORSE_descexpr4->nb, MORSE_descexpr4->bsiz,
			m, MORSE_descexpr4->n, 0, 0, MORSE_descexpr4->m, MORSE_descexpr4->n, 1, 1,
			morse_getaddr_cm, morse_getblkldd_cm, NULL );
	MORSE_Sequence_Wait(msequence);
	double total_loop_time =0.0;
	double loop_time = 0.0;
	for(p=0; p<nZmiss; p++)
	{

		lmiss.x[0]            = data->lmiss.x[p];
		lmiss.y[0]            = data->lmiss.y[p];

		START_TIMING(matrix_gen3)
			MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_desck_t, &data->lobs, &lmiss, &data->lm, truth_theta, data->dm, data->kernel_fun,  msequence, mrequest);
		MORSE_Sequence_Wait(msequence);
		STOP_TIMING(matrix_gen3);

		START_TIMING(matrix_gen4)
			MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_desck_a, &data->lobs, &lmiss, &data->lm, estimated_theta, data->dm, data->kernel_fun,  msequence, mrequest);
		MORSE_Sequence_Wait(msequence);
		STOP_TIMING(matrix_gen4);


		START_TIMING(copy2);
		//(5)Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile(MorseUpperLower, MORSE_desck_t, MORSE_desck_ttmp);
		VERBOSE(" Done.\n");
		STOP_TIMING(copy2);

		//  double *C = (double *) malloc(2 * n * sizeof(double));
		//  MORSE_Tile_to_Lapack( MORSE_desck_t, C, nZobs);
		//   print_dmatrix("testC", 2, 16, C, 2);
		//exit(0);

		START_TIMING(copy3);
		//(6)Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile(MorseUpperLower, MORSE_desck_a, MORSE_desck_atmp);
		VERBOSE(" Done.\n");
		STOP_TIMING(copy3);


		START_TIMING(loop_time);
		START_TIMING(trsm1);
		//(7) Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-1, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_a->n);
		STOP_TIMING(trsm1);

		START_TIMING(trsm2);
		//(8) Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a);
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_a->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(trsm2);

		START_TIMING(trsm3);
		//(9) Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
		VERBOSE("(9)Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t);
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_t->m, MORSE_desck_t->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(trsm3);

		START_TIMING(trsm4);
		//(10) Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
		VERBOSE("(10)Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t);
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_t->m, MORSE_desck_t->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(trsm4);

		START_TIMING(gevv1);
		//(12) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("(11)Calculate dgemm MORSE_descexpr4 = MORSE_desck_a^T * MORSE_desck_t... (Prediction Stage)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_t, 0, MORSE_descexpr3);
		flops = flops + FLOPS_DGEMM(MORSE_desck_ttmp->m, MORSE_desck_t->n, MORSE_descexpr3->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(gevv1);

		START_TIMING(gevv2);
		//(11) Calculate dgemm value= MORSE_desck_t^T * MORSE_desck_a
		VERBOSE("(12)Calculate dgemm MORSE_descexpr1 = MORSE_desck_t^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_a, 0, MORSE_descexpr1);
		flops = flops + FLOPS_DGEMM(MORSE_desck_ttmp->m, MORSE_desck_a->n, MORSE_descexpr1->n);	
		VERBOSE(" Done.\n");
		STOP_TIMING(gevv2);
		START_TIMING(gevv3);
		//(8) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_atmp
		VERBOSE("(13)Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_a, 0, MORSE_descexpr4);
		flops = flops + FLOPS_DGEMM(MORSE_desck_atmp->m, MORSE_desck_a->n, MORSE_descexpr4->n);	
		VERBOSE(" Done.\n");
		STOP_TIMING(gevv3);
		START_TIMING(gevv4);
		//(14) Calculate dgemm MORSE_desck_a= MORSE_descK_t * MORSE_desck_a (use k_t as k_a)
		VERBOSE("(14)Calculate dgemm MORSE_desck_a = MORSE_descK_ttmp * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile(MorseNoTrans, MorseNoTrans, 1, MORSE_descK_ttmp, MORSE_desck_a, 0, MORSE_desck_t);
		flops = flops + FLOPS_DGEMM(MORSE_desck_ttmp->m, MORSE_desck_a->n, MORSE_desck_t->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(gevv4);

		START_TIMING(trsm5);
		//(15) Triangular Solve (TRSM) k_atmp = TRSM(K_a^-1, k_atmp)
		VERBOSE("(15)Solving the linear system k_atmp = TRSM(K_t^-1, k_atmp) ...\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t);
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_t->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(trsm5);

		START_TIMING(trsm6);
		//(16) Triangular Solve (TRSM) k_a = TRSM(K_a^-T, k_a)
		VERBOSE("(16)Solving the linear system k_atmp = TRSM(K_t^-T, k_atmp) ...\n");
		MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t);
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_t->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(trsm6);
		START_TIMING(gevv5);
		//(13) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("(17)Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_t, 0, MORSE_descexpr2);
		flops = flops + FLOPS_DGEMM(MORSE_desck_atmp->m, MORSE_desck_t->n, MORSE_descexpr2->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(gevv5);


		STOP_TIMING(loop_time);

		total_loop_time +=loop_time;
		//stop
		//exit(0);
		//		double expr1=0;
		//		double expr2=0;
		//		double expr3=0;
		//		double expr4=0;


		//	MORSE_Tile_to_Lapack( MORSE_descexpr1, &expr1, 1);
		//	MORSE_Tile_to_Lapack( MORSE_descexpr2, &expr2, 1);
		//	MORSE_Tile_to_Lapack( MORSE_descexpr3, &expr3, 1);
		//	MORSE_Tile_to_Lapack( MORSE_descexpr4, &expr4, 1);

		MORSE_dgeadd_Tile(MorseNoTrans, 1 , MORSE_desctruthalpha, -2, MORSE_descexpr1);
		MORSE_dgeadd_Tile(MorseNoTrans, 1 , MORSE_descexpr1, 1, MORSE_descexpr2);

		MORSE_dgeadd_Tile(MorseNoTrans, 1 , MORSE_desctruthalpha,-1, MORSE_descexpr3);
		MORSE_dgeadd_Tile(MorseNoTrans, 1 , MORSE_descestimatedalpha, -1, MORSE_descexpr4);
		//		temp1 = truth_theta[0]-  2* data->expr1 + data->expr2;
		//		temp2 = truth_theta[0]-  data->expr3;
		//		temp3 = estimated_theta[0]- data->expr4;

		//		printf ("%f, %f, %f,%f, %f\n", matrix_gen1, matrix_gen2, copy1, cholesky1, cholesky2);
		//		printf ("%f, %f, %f, %f\n", matrix_gen3, matrix_gen4, copy2, copy3);
		//		printf ("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", trsm1, trsm2, trsm3, trsm3, trsm5, trsm6, gevv1, gevv2, gevv3, gevv4, gevv5);
		//*******************************

		//MORSE_Tile_to_Lapack( MORSE_descexpr2, temp1, m);
		//MORSE_Tile_to_Lapack( MORSE_descexpr3, temp2, m);
		//MORSE_Tile_to_Lapack( MORSE_descexpr4, temp3, m);	
		MORSE_dlacpy_Tile( MorseUpperLower, MORSE_descexpr2, B2);
		MORSE_dlacpy_Tile( MorseUpperLower, MORSE_descexpr3, B3);
		MORSE_dlacpy_Tile( MorseUpperLower, MORSE_descexpr4, B4);

		if(m==1)
		{
			printf("%f, %f, %f\n", temp1[0], temp2[0], temp3[0]);
			loe[p]  = temp1[0]/temp2[0]-1.0;
			mom[p]  = temp3[0]/temp1[0]-1.0;

		}
		if(m==2)
		{
			printf("%f, %f, %f, %f, %f, %f\n", temp1[0], temp1[2], temp2[0],temp2[2],temp3[0],temp3[2]);	
			loe[p]  = (temp1[0]+temp1[2])/(temp2[0]+temp2[2])-1.0;
			mom[p]  = (temp3[0]+temp3[2])/(temp1[0]+temp1[2])-1.0;
		}
		loe_sum += loe[p];
		mom_sum += mom[p];
	}

#if defined(CHAMELEON_USE_MPI)
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		printf("\nnZmiss = %d\n", nZmiss);
		printf("\nMLOE = %f", loe_sum/nZmiss);
		printf("\nMMOM = %f\n\n\n", mom_sum/nZmiss);
		data->mloe=(loe_sum/nZmiss);
		data->mmom=(mom_sum/nZmiss);
		printf(" ----(MLOE-MMOM) Gflop/s: %6.2f\n", flops / 1e9 / (total_loop_time  + cholesky1 +cholesky2));



#if defined(CHAMELEON_USE_MPI)
	}
#endif

	data->mloe=(loe_sum/nZmiss);
	data->mmom=(mom_sum/nZmiss);
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
	if ( MORSE_My_Mpi_Rank() == 0 )
	{
#endif
		//output
		results.mloe = data->mloe;
		results.mmom = data->mmom;
		results.mloe_exec = "sync";
		results.total_mloe_mmom_time = all_time;
		results.matrix_gen_mloe_mmom_time = matrix_gen1+ matrix_gen2;
		results.cho_fact_mloe_mmom_time = cholesky1 +cholesky2;
		results.loop_mloe_mmom_time = total_loop_time;
		results.total_mloe_mmom_flops = flops / 1e9 / (total_loop_time  + cholesky1 +cholesky2);
#if defined(CHAMELEON_USE_MPI)
	}
#endif

	fprintf(stderr," ---- mloe_mmom Time: %6.2f seconds\n\n", all_time);
}



void MORSE_dmle_mloe_mmom_Tile_Async(MLE_data *MORSE_data, double * truth_theta, double* estimated_theta, int nZmiss, int nZobs, int n)
	//! //Predict missing values base on a set of given values and covariance matrix
	/*!  -- MORSE-sync
	 * Returns the prediction Mean Square Error (MSE) as double
	 * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
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
	//truth_theta[0]=1; truth_theta[1]=0.1; truth_theta[2]=0.5;
	//estimated_theta[0]=1.01; estimated_theta[1]=0.09; estimated_theta[2]=0.49;
	printf("%f, %f, %f\n", truth_theta[0], truth_theta[1], truth_theta[2]);
	printf("%f, %f, %f\n", estimated_theta[0], estimated_theta[1], estimated_theta[2]);
	//exit(0);
	double loe_sum  = 0.0;
	double mom_sum  = 0.0;
	int i           = 0;
	int p           = 0;
	int v           = 0;
	int j           = 0;
	double all_time = 0.0;
	double cholesky1 = 0.0;
	double cholesky2 = 0.0;
	double matrix_gen1 = 0.0;
	double matrix_gen2 = 0.0;
	double matrix_gen3 = 0.0;
	double matrix_gen4 = 0.0;
	double copy1 = 0.0;
	double copy2 = 0.0;
	double copy3 = 0.0;

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
	double *loe=(double *) malloc(nZmiss * sizeof(double));
	double *mom=(double *) malloc(nZmiss * sizeof(double));

	MLE_data     *data                        = (MLE_data*) MORSE_data;
	MORSE_desc_t *MORSE_desck_t               = data->desck_t;
	MORSE_desc_t *MORSE_desck_a               = data->desck_a;
	MORSE_desc_t *MORSE_descK_t               = data->descK_t;
	MORSE_desc_t *MORSE_descK_a               = data->descK_a;
	MORSE_desc_t *MORSE_descK_ttmp            = data->descK_ttmp;
	MORSE_desc_t *MORSE_desck_atmp            = data->desck_atmp;
	MORSE_desc_t *MORSE_desck_ttmp            = data->desck_ttmp;
	MORSE_desc_t *MORSE_descexpr1             = data->descexpr1;
	MORSE_desc_t *MORSE_descexpr2             = data->descexpr2;
	MORSE_desc_t *MORSE_descexpr3             = data->descexpr3;
	MORSE_desc_t *MORSE_descexpr4             = data->descexpr4;
	MORSE_desc_t *MORSE_descestimatedalpha    = data->descestimatedalpha;
	MORSE_desc_t *MORSE_desctruthalpha        = data->desctruthalpha;
	MORSE_sequence_t *msequence               = (MORSE_sequence_t *)(data->sequence);
	MORSE_request_t *mrequest                 = (MORSE_request_t *)data->request;
	location lmiss;
	lmiss.x            = (double *) malloc(sizeof(double));
	lmiss.y            = (double *) malloc(sizeof(double));

	double* univariate_theta;
	double* univariate2_theta;
	double* univariate3_theta;
	double nu12;
	double rho;
	double sigma_square12;

	double flops = 0.0;
	START_TIMING(all_time);

	int m = MORSE_descestimatedalpha->m;
	double *truthalpha    = (double *) malloc( m * m * sizeof(double));
	double *estimatedalpha = (double *) malloc( m * m * sizeof(double));
	double *temp1 = (double *) malloc( m * m* sizeof(double));
	double *temp2 = (double *) malloc( m * m*sizeof(double));
	double *temp3 = (double *) malloc( m * m*sizeof(double));
	if(m ==1)
	{
		truthalpha[0]   = truth_theta[0];
		estimatedalpha[0] = estimated_theta[0];
	}

	if(m ==2)
	{
		double truth_nu12 = 0.5 * (truth_theta[3] + truth_theta[4]);
		double truth_rho = truth_theta[5] * sqrt( (tgamma(truth_theta[3] + 1)*tgamma(truth_theta[4] + 1)) /
				(tgamma(truth_theta[3]) * tgamma(truth_theta[4])) ) *
			tgamma(truth_nu12) / tgamma(truth_nu12 + 1);

		double estimated_nu12 = 0.5 * (estimated_theta[3] + estimated_theta[4]);
		double estimated_rho = estimated_theta[5] * sqrt( (tgamma(estimated_theta[3] + 1)*tgamma(estimated_theta[4] + 1)) /
				(tgamma(estimated_theta[3]) * tgamma(estimated_theta[4])) ) *
			tgamma(estimated_nu12) / tgamma(estimated_nu12 + 1);

		truthalpha[0]    = truth_theta[0];
		estimatedalpha[0] = estimated_theta[0];

		truthalpha[1]    = truthalpha[3] = truth_rho
			* sqrt(truth_theta[0] * truth_theta[1]);

		estimatedalpha[1] = estimatedalpha[3] = estimated_rho
			* sqrt(estimated_theta[0] * estimated_theta[1]);
		truthalpha[2]    = truth_theta[1];
		estimatedalpha[2] = estimated_theta[1];

	}
	MORSE_Lapack_to_Tile( truthalpha, m, MORSE_desctruthalpha);
	MORSE_Lapack_to_Tile( estimatedalpha, m, MORSE_descestimatedalpha);



	char * name= malloc(strlen(data->kernel_fun) + 1); 
	strcpy(name, data->kernel_fun);

	if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
	{
		printf("TODO:  running with Z=1 no more\n");
		data->kernel_fun= "bivariate_matern_parsimonious";
	}

	START_TIMING(matrix_gen1);
	VERBOSE("Create MORSE_descK_a Covariance Matrix (MLOE-MMOM).....");
	if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
	{

		univariate_theta =(double *) malloc(3 * sizeof(double));
		univariate2_theta =(double *) malloc(3 * sizeof(double));
		univariate3_theta =(double *) malloc(3 * sizeof(double));
		univariate_theta[0]=estimated_theta[0];
		univariate_theta[1]=estimated_theta[2];
		univariate_theta[2]=estimated_theta[3];

		MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, morse_desc_submatrix(MORSE_descK_a, 0,   0, MORSE_descK_a->m/2, MORSE_descK_a->n/2), &data->lobs, &data->lobs, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);


		nu12 = 0.5 * (estimated_theta[3] + estimated_theta[4]);
		rho = estimated_theta[5] * sqrt( (tgamma(estimated_theta[3] + 1)*tgamma(estimated_theta[4] + 1)) /
				(tgamma(estimated_theta[3]) * tgamma(estimated_theta[4])) ) *
			tgamma(nu12) / tgamma(nu12 + 1);
		sigma_square12 = rho * sqrt(estimated_theta[0]*estimated_theta[1]) ;

		univariate2_theta[0]=sigma_square12;
		univariate2_theta[1]=estimated_theta[2];
		univariate2_theta[2]=nu12;

		MORSE_MLE_dcmg_Tile_Async(MorseUpperLower,  morse_desc_submatrix(MORSE_descK_a, MORSE_descK_a->m/2,   0, MORSE_descK_a->m/2, MORSE_descK_a->n/2), &data->lobs, &data->lobs, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);
		//      MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descsubC21, &data->l1, &data->l1, &data->lm, univariate2_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

		univariate3_theta[0]=estimated_theta[1];
		univariate3_theta[1]=estimated_theta[2];
		univariate3_theta[2]=estimated_theta[4];

		MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, morse_desc_submatrix(MORSE_descK_a, MORSE_descK_a->m/2,   MORSE_descK_a->n/2, MORSE_descK_a->m/2, MORSE_descK_a->n/2), &data->lobs, &data->lobs, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);
	}
	else
	{               //(1)Generate the co-variance matrix descK_a
		MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descK_a,  &data->lobs, &data->lobs, &data->lm, estimated_theta, data->dm, data->kernel_fun,  msequence, mrequest);
	}
	VERBOSE(" Done.\n");
	MORSE_Sequence_Wait(msequence);
	STOP_TIMING(matrix_gen1);

	START_TIMING(matrix_gen2);
	//(2)Generate the co-variance matrix descK_t
	VERBOSE("Create MORSE_descK_t Covariance Matrix (MLOE-MMOM).....");
	MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_descK_t,  &data->lobs, &data->lobs, &data->lm, truth_theta, data->dm, data->kernel_fun,  msequence, mrequest);
	MORSE_Sequence_Wait(msequence);
	VERBOSE(" Done.\n");
	STOP_TIMING(matrix_gen2);


	START_TIMING(cholesky1);
	//(3)Cholesky factorization for the Co-variance matrix MORSE_descK_a
	VERBOSE("Cholesky factorization of MORSE_descK_a (MLOE-MMOM) .....");
	int success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descK_a, msequence, mrequest);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");
	STOP_TIMING(cholesky1);
	flops = flops + FLOPS_DPOTRF(MORSE_descK_a->m);

	START_TIMING(copy1);
	//(4)Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM)
	VERBOSE("Copy MORSE_descK_t Covariance Matrix to MORSE_descK_ttmp  (MLOE-MMOM).....");
	MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descK_t, MORSE_descK_ttmp, msequence, mrequest);
	VERBOSE(" Done.\n");
	STOP_TIMING(copy1);

	START_TIMING(cholesky2);
	flops = flops + FLOPS_DPOTRF(MORSE_descK_t->m);
	//(5)Cholesky factorization for the Co-variance matrix MORSE_descK_t
	VERBOSE("Cholesky factorization of MORSE_descK_t (MLOE-MMOM) .....");
	success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descK_t, msequence, mrequest);
	SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
	VERBOSE(" Done.\n");


	STOP_TIMING(cholesky2);
	MORSE_desc_t *B2;
	MORSE_desc_t *B3;
	MORSE_desc_t *B4;
	MORSE_Desc_Create_User( &B2, temp1, MorseComplexDouble, MORSE_descexpr2->mb, MORSE_descexpr2->nb, MORSE_descexpr2->bsiz,
			m, MORSE_descexpr2->n, 0, 0, MORSE_descexpr2->m, MORSE_descexpr2->n, 1, 1,
			morse_getaddr_cm, morse_getblkldd_cm, NULL );
	MORSE_Desc_Create_User( &B3, temp2, MorseComplexDouble, MORSE_descexpr3->mb, MORSE_descexpr3->nb, MORSE_descexpr3->bsiz,
			m, MORSE_descexpr3->n, 0, 0, MORSE_descexpr3->m, MORSE_descexpr3->n, 1, 1,
			morse_getaddr_cm, morse_getblkldd_cm, NULL );
	MORSE_Desc_Create_User( &B4, temp3, MorseComplexDouble, MORSE_descexpr4->mb, MORSE_descexpr4->nb, MORSE_descexpr4->bsiz,
			m, MORSE_descexpr4->n, 0, 0, MORSE_descexpr4->m, MORSE_descexpr4->n, 1, 1,
			morse_getaddr_cm, morse_getblkldd_cm, NULL );
	double total_loop_time =0.0;
	double loop_time = 0.0;
	for(p=0; p<nZmiss; p++)
	{

		//printf("p:%d\n", p);
		lmiss.x[0]            = data->lmiss.x[p];
		lmiss.y[0]            = data->lmiss.y[p];
		//printf("p:%d\n", p);
		START_TIMING(matrix_gen3)
			MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_desck_t, &data->lobs, &lmiss, &data->lm, truth_theta, data->dm, data->kernel_fun,  msequence, mrequest);
		//	MORSE_Sequence_Wait(msequence);
		STOP_TIMING(matrix_gen3);

		START_TIMING(matrix_gen4)
			MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, MORSE_desck_a, &data->lobs, &lmiss, &data->lm, estimated_theta, data->dm, data->kernel_fun,  msequence, mrequest);
		//	MORSE_Sequence_Wait(msequence);
		STOP_TIMING(matrix_gen4);

		//printf("p:%d\n", p);
		START_TIMING(copy2);
		//(5)Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_a to MORSE_descK_atmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_desck_t, MORSE_desck_ttmp, msequence, mrequest);
		VERBOSE(" Done.\n");
		STOP_TIMING(copy2);

		//  double *C = (double *) malloc(2 * n * sizeof(double));
		//  MORSE_Tile_to_Lapack( MORSE_desck_t, C, nZobs);
		//   print_dmatrix("testC", 2, 16, C, 2);
		//exit(0);

		START_TIMING(copy3);
		//(6)Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM)
		VERBOSE("Copy MORSE_desck_t to MORSE_desck_ttmp  (MLOE-MMOM).....");
		MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_desck_a, MORSE_desck_atmp, msequence, mrequest);
		VERBOSE(" Done.\n");
		STOP_TIMING(copy3);


		START_TIMING(loop_time);
		START_TIMING(trsm1);
		//(7) Triangular Solve (TRSM) k_a = TRSM(L_a^-1, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-1, k_a) ...(MLOE-MMOM)\n");


		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_a->n);
		STOP_TIMING(trsm1);

		START_TIMING(trsm2);
		//(8) Triangular Solve (TRSM) k_a = TRSM(L_a^-T, k_a)
		VERBOSE("Solving the linear system k_a = TRSM(L_a^-T, k_a) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_a, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_a->n);
		STOP_TIMING(trsm2);

		START_TIMING(trsm3);
		//(9) Triangular Solve (TRSM) k_t = TRSM(L_t^-1, k_t)
		VERBOSE("(9)Solving the linear system k_t = TRSM(L_t^-1, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_t->m, MORSE_desck_t->n);
		STOP_TIMING(trsm3);

		START_TIMING(trsm4);
		//(10) Triangular Solve (TRSM) k_t = TRSM(L_t^-T, k_t)
		VERBOSE("(10)Solving the linear system k_t = TRSM(L_a^-T, k_t) ...(MLOE-MMOM)\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_t, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_t->m, MORSE_desck_t->n);
		STOP_TIMING(trsm4);

		START_TIMING(gevv1);
		//(12) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("(11)Calculate dgemm MORSE_descexpr4 = MORSE_desck_a^T * MORSE_desck_t... (Prediction Stage)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_t, 0, MORSE_descexpr3, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DGEMM(MORSE_desck_ttmp->m, MORSE_desck_t->n, MORSE_descexpr3->n);
		STOP_TIMING(gevv1);

		START_TIMING(gevv2);
		//(11) Calculate dgemm value= MORSE_desck_t^T * MORSE_desck_a
		VERBOSE("(12)Calculate dgemm MORSE_descexpr1 = MORSE_desck_t^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_ttmp, MORSE_desck_a, 0, MORSE_descexpr1, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DGEMM(MORSE_desck_ttmp->m, MORSE_desck_a->n, MORSE_descexpr1->n);
		STOP_TIMING(gevv2);
		START_TIMING(gevv3);
		//(8) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_atmp
		VERBOSE("(13)Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (MLOE-MMOM)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_a, 0, MORSE_descexpr4, msequence, mrequest);

		flops = flops + FLOPS_DGEMM(MORSE_desck_atmp->m, MORSE_desck_a->n, MORSE_descexpr4->n);
		VERBOSE(" Done.\n");
		STOP_TIMING(gevv3);
		START_TIMING(gevv4);
		//(14) Calculate dgemm MORSE_desck_a= MORSE_descK_t * MORSE_desck_a (use k_t as k_a)
		VERBOSE("(14)Calculate dgemm MORSE_desck_a = MORSE_descK_ttmp * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile_Async(MorseNoTrans, MorseNoTrans, 1, MORSE_descK_ttmp, MORSE_desck_a, 0, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DGEMM(MORSE_desck_ttmp->m, MORSE_desck_a->n, MORSE_desck_t->n);
		STOP_TIMING(gevv4);

		START_TIMING(trsm5);
		//(15) Triangular Solve (TRSM) k_atmp = TRSM(K_a^-1, k_atmp)
		VERBOSE("(15)Solving the linear system k_atmp = TRSM(K_t^-1, k_atmp) ...\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_t->n);
		STOP_TIMING(trsm5);

		START_TIMING(trsm6);
		//(16) Triangular Solve (TRSM) k_a = TRSM(K_a^-T, k_a)
		VERBOSE("(16)Solving the linear system k_atmp = TRSM(K_t^-T, k_atmp) ...\n");
		MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descK_a, MORSE_desck_t, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DTRSM(MorseLeft, MORSE_descK_a->m, MORSE_desck_t->n);	
		STOP_TIMING(trsm6);
		START_TIMING(gevv5);
		//(13) Calculate dgemm value= MORSE_desck_a^T * MORSE_desck_t
		VERBOSE("(17)Calculate dgemm MORSE_descexpr1 = MORSE_desck_a^T * MORSE_desck_a... (Prediction Stage)");
		MORSE_dgemm_Tile_Async(MorseTrans, MorseNoTrans, 1, MORSE_desck_atmp, MORSE_desck_t, 0, MORSE_descexpr2, msequence, mrequest);
		VERBOSE(" Done.\n");
		flops = flops + FLOPS_DGEMM(MORSE_desck_atmp->m, MORSE_desck_t->n, MORSE_descexpr2->n);
		STOP_TIMING(gevv5);

		//stop
		//exit(0);
		//              double expr1=0;
		//              double expr2=0;
		//              double expr3=0;
		//              double expr4=0;


		//      MORSE_Tile_to_Lapack( MORSE_descexpr1, &expr1, 1);
		//      MORSE_Tile_to_Lapack( MORSE_descexpr2, &expr2, 1);
		//      MORSE_Tile_to_Lapack( MORSE_descexpr3, &expr3, 1);
		//      MORSE_Tile_to_Lapack( MORSE_descexpr4, &expr4, 1);



		STOP_TIMING(loop_time);

		total_loop_time +=loop_time;

		MORSE_dgeadd_Tile_Async(MorseNoTrans, 1 , MORSE_desctruthalpha, -2, MORSE_descexpr1, msequence, mrequest);
		MORSE_dgeadd_Tile_Async(MorseNoTrans, 1 , MORSE_descexpr1, 1, MORSE_descexpr2, msequence, mrequest);

		MORSE_dgeadd_Tile_Async(MorseNoTrans, 1 , MORSE_desctruthalpha,-1, MORSE_descexpr3, msequence, mrequest);
		MORSE_dgeadd_Tile_Async(MorseNoTrans, 1 , MORSE_descestimatedalpha, -1, MORSE_descexpr4, msequence, mrequest);
		//              temp1 = truth_theta[0]-  2* data->expr1 + data->expr2;
		//              temp2 = truth_theta[0]-  data->expr3;
		//              temp3 = estimated_theta[0]- data->expr4;

		//              printf ("%f, %f, %f,%f, %f\n", matrix_gen1, matrix_gen2, copy1, cholesky1, cholesky2);
		//              printf ("%f, %f, %f, %f\n", matrix_gen3, matrix_gen4, copy2, copy3);
		//              printf ("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", trsm1, trsm2, trsm3, trsm3, trsm5, trsm6, gevv1, gevv2, gevv3, gevv4, gevv5);
		//*******************************

		//		MORSE_Tile_to_Lapack( MORSE_descexpr2, temp1, m);
		//		MORSE_Tile_to_Lapack( MORSE_descexpr3, temp2, m);
		//		MORSE_Tile_to_Lapack( MORSE_descexpr4, temp3, m);

		MORSE_dlacpy_Tile_Async( MorseUpperLower, MORSE_descexpr2, B2, msequence, mrequest);
		MORSE_dlacpy_Tile_Async( MorseUpperLower, MORSE_descexpr3, B3, msequence, mrequest);
		MORSE_dlacpy_Tile_Async( MorseUpperLower, MORSE_descexpr4, B4, msequence, mrequest);

		MORSE_Sequence_Wait(msequence);
		if(m==1)
		{ 
			printf("%f, %f, %f, \n", temp1[0], temp2[0], temp3[0]);
			loe_sum +=  temp1[0]/temp2[0]-1.0;
			mom_sum += temp3[0]/temp1[0]-1.0;

		}
		if(m==2)
		{
			// loe=(double *) malloc(nZmiss * sizeof(double));
			// mom=(double *) malloc(nZmiss * sizeof(double));
			printf("%f, %f, %f,\n", temp1[0],  temp2[0], temp3[0]);
			loe_sum += (temp1[0]+temp1[2])/(temp2[0]+temp2[2])-1.0;
			mom_sum += (temp3[0]+temp3[2])/(temp1[0]+temp1[2])-1.0;
			//printf("%f, %f, %f, %f, %f, %f\n", temp1[0], temp1[2], temp2[0],temp2[2],temp3[0],temp3[2]);
		}
		//loe_sum += loe[p];
		//mom_sum += mom[p];
		//
		//	printf("%f, %f, %f, %f, %f, %f\n", temp1[0], temp1[2], temp2[0],temp2[2],temp3[0],temp3[2]);
	}
#if defined(CHAMELEON_USE_MPI)
	if(MORSE_My_Mpi_Rank() == 0)
	{
#endif
		printf("\nnZmiss = %d\n", nZmiss);
		printf("\nMLOE = %f", loe_sum/nZmiss);
		printf("\nMMOM = %f\n\n\n", mom_sum/nZmiss);
		data->mloe=(loe_sum/nZmiss);
		data->mmom=(mom_sum/nZmiss);
		printf(" ----(MLOE-MMOM) Gflop/s: %6.2f\n", flops / 1e9 / (total_loop_time  + cholesky1 +cholesky2));
#if defined(CHAMELEON_USE_MPI)
	}
#endif
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
	if ( MORSE_My_Mpi_Rank() == 0 )
	{
#endif
		//output
		results.mloe = data->mloe;
		results.mmom = data->mmom;
		results.mloe_exec = "sync";
		results.total_mloe_mmom_time = all_time;
		results.matrix_gen_mloe_mmom_time = matrix_gen1+ matrix_gen2;
		results.cho_fact_mloe_mmom_time = cholesky1 +cholesky2;
		results.loop_mloe_mmom_time = total_loop_time;
		results.total_mloe_mmom_flops = flops / 1e9 / (total_loop_time  + cholesky1 +cholesky2);
#if defined(CHAMELEON_USE_MPI)
	}
#endif

	fprintf(stderr," ---- mloe_mmom Time: %6.2f seconds\n\n", all_time);
}








//init Chameleon descriptors
void MORSE_dmle_Call(MLE_data  *data, int ncores,int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
	//! //Initiate MORSE and allocate different descriptors for
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

	MORSE_sequence_t *msequence;
	MORSE_request_t mrequest[2] = { MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER };
	MORSE_desc_t *MORSE_descC	 = NULL;
	MORSE_desc_t *MORSE_descsubC11   = NULL;
	MORSE_desc_t *MORSE_descsubC12   = NULL;
	//	MORSE_desc_t *MORSE_descsubC21   = NULL;
	MORSE_desc_t *MORSE_descsubC22   = NULL;
	MORSE_desc_t *MORSE_descZ	 = NULL;
	MORSE_desc_t *MORSE_descZ1       = NULL;
	MORSE_desc_t *MORSE_descZ2       = NULL;
	MORSE_desc_t *MORSE_descZcpy	 = NULL;
	MORSE_desc_t *MORSE_descproduct	 = NULL;
	MORSE_desc_t *MORSE_descproduct1 = NULL;
	MORSE_desc_t *MORSE_descproduct2 = NULL;
	MORSE_desc_t *MORSE_descdet	 = NULL;
	//MORSE_desc_t *MORSE_descZmiss	 = NULL;
	//MORSE_desc_t *MORSE_descC12	 = NULL;
	//MORSE_desc_t *MORSE_descC22	 = NULL;
	//MORSE_desc_t *MORSE_descmse	 = NULL;
	//MORSE_desc_t *MORSE_descZactual	= NULL;
	//MORSE_desc_t *MORSE_descZobs	 = NULL;



	// For ditributed system and should be removed
	double *Zcpy = (double *) malloc(N * sizeof(double));

	//Identifies a set of routines sharing common exception handling.
	MORSE_Sequence_Create(&msequence);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC,NULL , MorseRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ, NULL, MorseRealDouble, dts, dts, dts * dts, N, 1,  0, 0, N, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZcpy, Zcpy, MorseRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descproduct, &data->dotp, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descproduct1, &data->dotp1, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descproduct2, &data->dotp2, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descdet, &data->det, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);


	//printf("=====%d\n", MORSE_descZ->m/2);

	//Two submatrix descriptors for bi-variate profiling 
	//	MORSE_descZ1 = morse_desc_submatrix(MORSE_descZ,   0,   0, MORSE_descZ->m/2, 1   );
	//	MORSE_descZ2 = morse_desc_submatrix(MORSE_descZ,   MORSE_descZ->m/2,   0, MORSE_descZ->m/2,   1);

	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ1, NULL, MorseRealDouble, dts, dts, dts * dts, N/2, 1, 0, 0, N/2, 1, p_grid, q_grid);
	EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ2, NULL, MorseRealDouble, dts, dts, dts * dts, N/2, 1, 0, 0, N/2, 1, p_grid, q_grid);


	MORSE_descsubC11 = morse_desc_submatrix(MORSE_descC, 0,   0, MORSE_descC->m/2, MORSE_descC->n/2);
	MORSE_descsubC12 = morse_desc_submatrix(MORSE_descC, MORSE_descC->m/2,   0, MORSE_descC->m/2, MORSE_descC->n/2);
	//MORSE_descsubC21 = morse_desc_submatrix(MORSE_descC, 0,   MORSE_descC->m/2, MORSE_descC->m/2, MORSE_descC->n/2);
	MORSE_descsubC22 = morse_desc_submatrix(MORSE_descC, MORSE_descC->m/2,   MORSE_descC->n/2, MORSE_descC->m/2, MORSE_descC->n/2);

	//if(nZmiss != 0)
	//{
	//      if(strcmp(data->actualZFPath,"")==0)
	//	{
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, &Zcpy[nZmiss], MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, Zcpy, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);	
	//	}
	//	else
	//	{
	//		MORSE_descZobs = MORSE_descZcpy;
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
	//	}

	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);
	//		EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
	//	}
	//Fill data struct
	data->descC		= MORSE_descC;
	data->descsubC11        = MORSE_descsubC11;
	data->descsubC12        = MORSE_descsubC12;
	//	data->descsubC21        = MORSE_descsubC21;
	data->descsubC22        = MORSE_descsubC22;
	data->descZ		= MORSE_descZ;
	data->descZ1            = MORSE_descZ1;
	data->descZ2            = MORSE_descZ2;
	data->descZcpy		= MORSE_descZcpy;
	data->descdet		= MORSE_descdet;
	data->descproduct	= MORSE_descproduct;
	data->descproduct1      = MORSE_descproduct1;
	data->descproduct2      = MORSE_descproduct2;
	//data->descZmiss		= MORSE_descZmiss;
	//data->descC12		= MORSE_descC12;
	//data->descC22		= MORSE_descC22;
	//data->descmse		= MORSE_descmse;
	//data->descZactual	= MORSE_descZactual;
	//data->descZobs		= MORSE_descZobs;
	data->sequence		= msequence;
	data->request		= mrequest;
	//stop gsl error handler
	gsl_set_error_handler_off () ;

}


