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
#include "../include/diag.h"
#include "../include/MLE_approx.h"
//***************************************************************************************

void MORSE_MLE_dzvg_diag_Tile (MLE_data *data,  double * Nrand, double * initial_theta, int n, int ts, int log)
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- MORSE-sync 
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     * 	                     that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] ts: tile size (MB) is used only in the case of HiCMA not MORSE.
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

    //In the case of testing mode, Z should be generated using Nrand and initial_theta
    //if (test == 1)    
    //	{
    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase - diagonal approximation).....");
    //	MORSE_dlaset_Tile(MorseUpperLower, 0, 1, data->descC);
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);


    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
    {

        univariate_theta =(double *) malloc(3 * sizeof(double));
        univariate2_theta =(double *) malloc(3 * sizeof(double));
        univariate3_theta =(double *) malloc(3 * sizeof(double));
        univariate_theta[0]=initial_theta[0];
        univariate_theta[1]=initial_theta[2];
        univariate_theta[2]=initial_theta[3];

        MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, data->descsubC11, &data->l1, &data->l1, &data->lm, univariate_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);


        nu12 = 0.5 * (initial_theta[3] + initial_theta[4]);
        rho = initial_theta[5] * sqrt( (tgamma(initial_theta[3] + 1)*tgamma(initial_theta[4] + 1)) /
                (tgamma(initial_theta[3]) * tgamma(initial_theta[4])) ) *
            tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(initial_theta[0]*initial_theta[1]) ;

        univariate2_theta[0]=sigma_square12;
        univariate2_theta[1]=initial_theta[2];
        univariate2_theta[2]=nu12;

        MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, data->descsubC12, &data->l1, &data->l1, &data->lm, univariate2_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

        //      MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descsubC21, &data->l1, &data->l1, &data->lm, univariate2_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);

        univariate3_theta[0]=initial_theta[1];
        univariate3_theta[1]=initial_theta[2];
        univariate3_theta[2]=initial_theta[4];
        MORSE_MLE_dcmg_Tile_Async(MorseUpperLower, data->descsubC22, &data->l1, &data->l1, &data->lm, univariate3_theta, data->dm, "univariate_matern_stationary",  msequence, &mrequest[0]);
    }
    else    if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
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
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase - diagonal approximation) .....");
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase- diagonal approximation) .....");
    int success = MORSE_dpotrf_Tile(MorseLower, data->descC);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication    
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase - diagonal approximation) .....");
    MORSE_dtrmm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if(log==1)
    {
        double *z;
        MORSE_desc_t *MORSE_descZ = (MORSE_desc_t *)(data->descZ);
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....");
#if defined(CHAMELEON_USE_MPI)
        z = (double *) malloc(n * sizeof(double));
        MORSE_Tile_to_Lapack( MORSE_descZ, z, n);
        if ( MORSE_My_Mpi_Rank() == 0 )
            write_vectors(z, data, n);
        free(z);
#else
        z = MORSE_descZ->mat;
        write_vectors(z, data, n);
        free(z);
#endif
        VERBOSE(" Done.\n");
    }
    //	}
    /*	else
        {
        double * streamdata;
        streamdata=(double *) malloc(n * sizeof(double));

    //Reading Observations from disk and copy it to data->descZ
    VERBOSE("Reading Observations from disk .....- diagonal approximation");
    streamdata = readObsFile(data->obsFPath, n);        
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(data->sequence);
    VERBOSE(" Done.\n");
    free(streamdata);
    }
    */
    MORSE_dlaset_Tile(MorseUpperLower, 0,0, data->descC);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous- diagonal approximation)\n");
    VERBOSE("************************************************************\n");
}


void MORSE_MLE_dzvg_diag_Tile_Async(MLE_data *data,  double * Nrand, double * initial_theta, int n, int ts, int log)
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- MORSE-Async
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] ts: tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * */
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    //In the case of testing mode, Z should be generated using Nrand and initial_theta
    //if (test ==1)
    //{
    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    //MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, data->descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm, data->kernel_fun,  msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
    int success = MORSE_dpotrf_Tile_Async(MorseLower, data->descC,  msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
    MORSE_dtrmm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ, msequence, mrequest);
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

    //}
    /*  else
        {
        double * streamdata;
        streamdata=(double *) malloc(n * sizeof(double));

    //Reading Observations from disk and copy it to data->descZ
    VERBOSE("Reading Observations from disk .....");
    streamdata = readObsFile(data->obsFPath, n);
    MORSE_MLE_dzcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    VERBOSE(" Done.\n");
    free(streamdata);
    }
    */
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous)\n");
    VERBOSE("************************************************************\n");
}



double MORSE_dmle_diag_Tile(unsigned n, const double * theta, double * grad, void * MORSE_data) {
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
    int N, NRHS, success, diag_thick;
    double flops =0.0;	
    double* univariate_theta;
    double* univariate2_theta;
    double* univariate3_theta;
    double nu12;
    double rho;
    double sigma_square12;
    int num_params = 0;
    MLE_data* data	= ((MLE_data*)MORSE_data);
    data->det	= 0;
    data->dotp	= 0;
    diag_thick = data->diag_thick;
    int i = 0;
    MORSE_desc_t *MORSE_descC        = (MORSE_desc_t *) data->descC;
    MORSE_desc_t *MORSE_descsubC11   = (MORSE_desc_t *) data->descsubC11;
    MORSE_desc_t *MORSE_descsubC12   = (MORSE_desc_t *) data->descsubC12;
    //MORSE_desc_t *MORSE_descsubC21   = (MORSE_desc_t *) data->descsubC21;
    MORSE_desc_t *MORSE_descsubC22   = (MORSE_desc_t *) data->descsubC22;
    MORSE_desc_t *MORSE_descZ        = (MORSE_desc_t *) data->descZ;
    MORSE_desc_t *MORSE_descZ1       = (MORSE_desc_t *) data->descZ1;
    MORSE_desc_t *MORSE_descZ2       = (MORSE_desc_t *) data->descZ2;
    MORSE_desc_t *MORSE_descZcpy     = (MORSE_desc_t *) data->descZcpy;
    MORSE_desc_t *MORSE_descdet      = (MORSE_desc_t *) data->descdet;
    MORSE_desc_t *MORSE_descproduct  = (MORSE_desc_t *) data->descproduct;
    MORSE_desc_t *MORSE_descproduct1 = (MORSE_desc_t *) data->descproduct1;
    MORSE_desc_t *MORSE_descproduct2 = (MORSE_desc_t *) data->descproduct2;
    MORSE_sequence_t *msequence      = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest       = (MORSE_request_t *) data->request;

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
    else
    {	
        VERBOSE("Re-store the original Z vector...");
        MORSE_dlacpy_Tile(MorseUpperLower ,MORSE_descZcpy,MORSE_descZ);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(dzcpy_time);	

    MORSE_dlaset_Tile(MorseUpperLower, 0, 0, MORSE_descC);
    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix...");
    START_TIMING(matrix_gen_time);

    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
    {
        printf("TODO: Generate New Covariance Matrix with bivariate_matern_parsimonious2\n");
        return;

    }


    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC, &data->l1, &data->l1, &data->lm, (double *)theta, data->dm, data->kernel_fun, msequence, &mrequest[0]);


    //	MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC, msequence, &mrequest[0], &data->l1, &data->l1, (double *)theta,  data->dm, diag_thick);    
    MORSE_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");


    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma...");
    START_TIMING(time_facto);
    success = MORSE_dpotrf_diag_Tile(MorseLower, MORSE_descC, diag_thick);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_DPOTRF(N);
    VERBOSE(" Done.\n");

    //       double *C = (double *) malloc(N * N * sizeof(double));
    //      MORSE_Tile_to_Lapack( data->descC, C, N);

    //double sum=0;
    //for(i=0;i<N*N;i++)
    //      sum+=C[i];
    //printf("sum= %f\n", sum);
    //exit(0);
    /*
       double *C = (double *) malloc(N * N * sizeof(double));
       MORSE_Tile_to_Lapack( MORSE_descC, C, N);
       print_dmatrix("test", 5, 5, C, 5);
       exit(0);
       */
    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant ...");
    START_TIMING(logdet_calculate);
    MORSE_MLE_dmdet_Tile_Async(MORSE_descC, msequence, &mrequest[0], MORSE_descdet);
    MORSE_Sequence_Wait(msequence);
    logdet= 2*data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system ...\n");
    START_TIMING(time_solve);
    MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_DTRSM(MorseLeft,N, NRHS);
    VERBOSE(" Done.\n");    


    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...");
    //dotp=0;
    //MORSE_MLE_core_ddotp_Async(MORSE_descZ,MORSE_descproduct,msequence, &mrequest[0]);
    //MORSE_Sequence_Wait(msequence);        
    MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct); 


    if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0 )
    {

        loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
        MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ1, MORSE_descZ1, 0, MORSE_descproduct1);
        MORSE_dgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ2, MORSE_descZ2, 0, MORSE_descproduct2);
        data->variance1 = (1.0/(N/2)) * data->dotp1;
        data->variance2 = (1.0/(N/2)) * data->dotp2;
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



    }
    else
    {

        loglik = -0.5 * data->dotp -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
        data->variance= theta[0];
    }
    VERBOSE(" Done.\n");

    //loglik = -0.5 * data->dotp -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
    // loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->dotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
    // VERBOSE(" Done.\n");

    //multiplicative scale
    // data->variance = (1.0/N) * data->dotp;
    //data->variance = theta[0];

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    // MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(MORSE_My_Mpi_Rank() == 0)
    {
#endif
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
        printf(" ---- Total Time: %6.2f\n", matrix_gen_time+ time_facto + logdet_calculate + time_solve);
        //fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );
        printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter+=/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter+=flops / 1e9 / (time_facto +time_solve);
    data->final_loglik=loglik;

    //output
    results.final_loglik=loglik;


    return loglik;
}

double MORSE_dmle_diag_Tile_Async(unsigned n, const double * theta, double * grad, void * MORSE_data) {
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
    int N, NRHS, success, diag_thick;

    MLE_data* data	= ((MLE_data*)MORSE_data);
    data->det	= 0;
    data->dotp	= 0;
    diag_thick = data->diag_thick;
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
    //	MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest, &data->l1, &data->l1,(double*) theta,  data->dm, diag_thick);    
    MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC, &data->l1, &data->l1, &data->lm, (double *)theta, data->dm, data->kernel_fun,  diag_thick, msequence, &mrequest[0]);


    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");
    MORSE_Sequence_Wait(msequence);

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma...");
    START_TIMING(time_facto);
    success = MORSE_dpotrf_diag_Tile_Async(MorseLower, MORSE_descC, diag_thick, msequence, mrequest);
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
        //fprintf(stderr,"------ddotproduct: %2.6f ", data->dotp);
        //fprintf(stderr,"------logdet: %2.6f ", logdet);
        //fprintf(stderr,"------det: %.*e ", det);
        //fprintf(stderr,"------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
        //fprintf(stderr," ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----LogLi: %2.6f\n", theta[0], theta[1], theta[2],loglik);
        //reformat
        fprintf(stderr," %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n", data->iter_count+1,  theta[0], theta[1], theta[2],loglik);

        if(data->log == 1)
            fprintf(data->pFileLog, " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n", data->iter_count+1,  theta[0], theta[1], theta[2],loglik);

        //fprintf(stderr," ---- Facto Time: %6.2f\n", time_facto);
        //fprintf(stderr," ---- logdet Time: %6.2f\n", logdet_calculate);
        //fprintf(stderr," ---- dtrsm Time: %6.2f\n", time_solve);
        //fprintf(stderr," ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
        //fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", dzcpy_time);
        //fprintf(stderr," ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
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
/*
   void MORSE_MLE_diag_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid, int mse_flag)
   {

   MORSE_desc_t *MORSE_descZmiss   = NULL;
   MORSE_desc_t *MORSE_descC12     = NULL;
   MORSE_desc_t *MORSE_descC22     = NULL;
   MORSE_desc_t *MORSE_descmse     = NULL;
   MORSE_desc_t *MORSE_descZactual = NULL;
   MORSE_desc_t *MORSE_descZobs    = NULL;
   MLE_data     *data              = (MLE_data*) MORSE_data;

   if(nZmiss <= 0)
   {
   fprintf(stderr," Number of missing values should be positive value\n");
   return -1;
   }
//Descriptors Creation
EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
if( mse_flag == 1)
{
EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);

EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
}
EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

//Initiate data descriptors
data->descZmiss         = MORSE_descZmiss;
data->descC12           = MORSE_descC12;
data->descC22           = MORSE_descC22;
data->descmse           = MORSE_descmse;
data->descZactual       = MORSE_descZactual;
data->descZobs          = MORSE_descZobs;
}
*/
double MORSE_dmle_diag_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n)
    //! //Predict missing values base on a set of given values and covariance matrix
    /*!  -- MORSE-sync
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
    //double *z = NULL,
    //double  *streamdata		= NULL;
    double time_solve 		= 0.0;
    double time_gemm 		= 0.0;
    double mat_gen_time 		= 0.0;
    double time_mse 		= 0.0;
    double flops 			= 0.0;

    MORSE_desc_t *MORSE_descZmiss   = NULL;
    MORSE_desc_t *MORSE_descC12     = NULL;
    MORSE_desc_t *MORSE_descC22     = NULL;
    MORSE_desc_t *MORSE_descmse     = NULL;
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

    //Descriptors Creation
    // EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
    //if( Zactual != NULL)
    //{
    //      EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
    //    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    // }
    // EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
    //EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);


    //Initiate data descriptors
    MORSE_descZmiss         = data->descZmiss;
    MORSE_descC12           = data->descC12;
    MORSE_descC22           = data->descC22;
    MORSE_descmse           = data->descmse;
    MORSE_descZactual       = data->descZactual;
    MORSE_descZobs          = data->descZobs;

    int diag_thick = data->diag_thick;


    MORSE_dlaset_Tile(MorseUpperLower, 0,0, MORSE_descC12);
    MORSE_dlaset_Tile(MorseUpperLower, 0,0, MORSE_descC22);

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

    //      for (i=0;i<100;i++)
    //      printf("%f, %f, %f\n", data->lobs.x[i], data->lobs.y[i], Zobs[i]);

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
    theta[0]=data->variance;
    printf("estimated parameters: %f - %f - %f\n", theta[0], theta[1], theta[2]);


    START_TIMING(mat_gen_time);
    //Generate C22 covariance matrix
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
    //MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  &data->lobs, &data->lobs, theta, data->dm, data->diag_thick);

    MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  data->diag_thick, msequence, mrequest);


    MORSE_Sequence_Wait(msequence); 
    //flops = flops + FLOPS_DPOTRF(nZobs);
    VERBOSE(" Done.\n");


    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
    //	MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  &data->lmiss, &data->lobs, theta, data->dm, data->diag_thick);

    MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC22, &data->lmiss, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  diag_thick, msequence, mrequest);

    //flops = flops + FLOPS_DPOTRF(nZmiss);
    MORSE_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);


    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
    //MORSE_dposv_Tile(MorseLower, MORSE_descC22, MORSE_descZobs);
    MORSE_dpotrf_diag_Tile(MorseLower, MORSE_descC22, diag_thick);
    MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC22, MORSE_descZobs);
    MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseTrans, MorseNonUnit, 1, MORSE_descC22, MORSE_descZobs);
    flops = flops + FLOPS_DPOTRF(nZobs);
    flops = flops + FLOPS_DTRSM(MorseLeft, nZobs, nZobs);
    STOP_TIMING(time_solve);


    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
    MORSE_dgemm_Tile (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss);
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //return back descZmiss to zmiss vector
    MORSE_Tile_to_Lapack( MORSE_descZmiss, Zmiss, nZmiss);

    //Estimate Mean Square Error
    if( Zactual != NULL)
    {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
        MORSE_MLE_dmse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
        VERBOSE(" Done.\n");
        MORSE_Sequence_Wait(msequence);
        STOP_TIMING(time_mse);
        data->mserror /= nZmiss;
    }
    else
        data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
    if(MORSE_My_Mpi_Rank() == 0)
    {
#endif
        if(data->log == 1)
            fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

        write_prediction_result("predict_result.dat", n, diag_thick, 0, 0, data->mserror, (mat_gen_time+time_solve+ time_gemm), (flops / 1e9 / (time_solve)));

#if defined(CHAMELEON_USE_MPI)
    }
#endif


    return data->mserror;


}



double MORSE_dmle_Predict_diag_Tile_Async(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, int n)
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
    double *streamdata 		= NULL;
    double mat_gen_time 		= 0.0;
    double time_solve 		= 0.0;
    double time_mse 		= 0.0;
    double flops 			= 0.0; 
    MLE_data * data			= (MLE_data*)MORSE_data;
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
        shuffle(z,  &data->l1, n);

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
        temp_loc.x=&l1->x[nZmiss];
        temp_loc.y=&l1->y[nZmiss];
        l2 = &temp_loc;

        //l1=(location *) calloc((size_t)nZmiss , sizeof(location));

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
    //	MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  l2, l2, theta, data->dm, data->diag_thick);
    MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC22, l2, l2, &data->lm, theta, data->dm, data->kernel_fun,  data->diag_thick, msequence, mrequest);

    //flops = flops + FLOPS_DPOTRF(nZobs);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
    //	MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  l1, l2, theta, data->dm, data->diag_thick);

    //	MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC12, l1, l2, &data->lm, theta, data->dm, data->kernel_fun,  data->diag_thick, msequence);

    MORSE_MLE_dcmg_diag_Tile_Async(MorseLower, MORSE_descC12, l1, l2, &data->lm, theta, data->dm, data->kernel_fun,  data->diag_thick, msequence, mrequest);


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
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
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
            fprintf(data->pFileLog, "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n", nZmiss,  (mat_gen_time+time_solve+time_mse), (flops / 1e9 / (time_solve)), data->mserror );

        // write_prediction_result("predict_result.dat", n, nZmiss, data->mserror, (mat_gen_time+time_solve+ time_mse), (flops / 1e9 / (time_solve )));

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;
}



//init Chameleon descriptors
void MORSE_dmle_diag_Call(MLE_data  *data, int ncores,int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
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
    MORSE_desc_t *MORSE_descC       = NULL;
    MORSE_desc_t *MORSE_descZ       = NULL;
    MORSE_desc_t *MORSE_descZcpy    = NULL;
    MORSE_desc_t *MORSE_descproduct = NULL;
    MORSE_desc_t *MORSE_descdet     = NULL;
    //MORSE_desc_t *MORSE_descZmiss = NULL;
    //MORSE_desc_t *MORSE_descC12   = NULL;
    //MORSE_desc_t *MORSE_descC22   = NULL;
    //MORSE_desc_t *MORSE_descmse   = NULL;
    //MORSE_desc_t *MORSE_descZactual       = NULL;
    //MORSE_desc_t *MORSE_descZobs  = NULL;



    // For ditributed system and should be removed
    double *Zcpy = (double *) malloc(N * sizeof(double));

    //Identifies a set of routines sharing common exception handling.
    MORSE_Sequence_Create(&msequence);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC,NULL , MorseRealDouble, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ, NULL, MorseRealDouble, dts, dts, dts * dts, N, 1,  0, 0, N, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZcpy, Zcpy, MorseRealDouble, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descproduct, &data->dotp, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descdet, &data->det, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    //if(nZmiss != 0)
    //{
    //      if(strcmp(data->actualZFPath,"")==0)
    //      {
    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, &Zcpy[nZmiss], MorseRealDouble, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, Zcpy, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
    //      }
    //      else
    //      {
    //              MORSE_descZobs = MORSE_descZcpy;
    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealDouble, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);
    //      }

    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealDouble, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealDouble, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);
    //              EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    //      }
    //Fill data struct
    data->descC             = MORSE_descC;
    data->descZ             = MORSE_descZ;
    data->descZcpy          = MORSE_descZcpy;
    data->descdet           = MORSE_descdet;
    data->descproduct       = MORSE_descproduct;
    //data->descZmiss       = MORSE_descZmiss;
    //data->descC12         = MORSE_descC12;
    //data->descC22         = MORSE_descC22;
    //data->descmse         = MORSE_descmse;
    //data->descZactual     = MORSE_descZactual;
    //data->descZobs                = MORSE_descZobs;
    data->sequence          = msequence;
    data->request           = mrequest;
    //stop gsl error handler
    gsl_set_error_handler_off () ;

}


