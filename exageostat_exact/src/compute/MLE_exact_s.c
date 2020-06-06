/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
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
 * @date 2020-01-19
 *
 **/
#include "../include/MLE_exact_s.h"
//***************************************************************************************
void MORSE_MLE_szvg_Tile (MLE_data *data,  float * Nrand, double * initial_theta, int n, int dts, int log)
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- MORSE-sync 
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     * 	                     that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] dts: tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] test: if 0 -> real data mode, 1 ->test data mode.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * */
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    //In the case of testing mode, Z should be generated using Nrand and initial_theta
    //if (test == 1)    
    //{
    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase - single precision) .....");
    //MORSE_MLE_scmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    MORSE_MLE_scmg_Tile_Async(MorseLower, data->descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm, data->kernel_fun,  msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase - single precision) .....");
    MORSE_MLE_szcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase - single precision) .....");
    int success = MORSE_spotrf_Tile(MorseLower, data->descC);
    //printf(" success=%d \n", success);
    //exit(0);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication    
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase - single precision) .....");
    MORSE_strmm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if(log==1)
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

    /*	}
        else
        {
        double * streamdata;
        streamdata=(double *) malloc(n * sizeof(double));

    //Reading Observations from disk and copy it to data->descZ
    VERBOSE("Reading Observations from disk .....");
    streamdata = readObsFile(data->obsFPath, n);        
    MORSE_MLE_szcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(data->sequence);
    VERBOSE(" Done.\n");
    free(streamdata);
    }
    */
    MORSE_slaset_Tile(MorseUpperLower, 0, 0, data->descC);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous - single precision)\n");
    VERBOSE("************************************************************\n");
}



void MORSE_MLE_szcpy( MLE_data *data, double *streamdata)
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    VERBOSE("Copy Z from vector to decriptor.\n");
    MORSE_MLE_szcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    VERBOSE("Done Z copying step.\n");
    VERBOSE("************************************************************\n");
}



void MORSE_MLE_szvg_Tile_Async(MLE_data *data,  float * Nrand, double * initial_theta, int n, int dts, int log)
    //! Generate Observations Vector (Z) for testing Maximum
    /*! Likelihood function -- MORSE-Async
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] dts: tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] test: if 0 -> real data mode, 1 ->test data mode.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * */
{
    MORSE_sequence_t *msequence     = (MORSE_sequence_t *) data->sequence;
    MORSE_request_t  *mrequest      = (MORSE_request_t *) data->request;
    //In the case of testing mode, Z should be generated using Nrand and initial_theta
    //       if (test ==1)
    //      {
    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase - single precision).....");
    //MORSE_MLE_scmg_Tile_Async(MorseLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    MORSE_MLE_scmg_Tile_Async(MorseLower, data->descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm, data->kernel_fun,  msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase - single precision) .....");
    MORSE_MLE_szcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase - single precision) .....");
    int success = MORSE_spotrf_Tile_Async(MorseLower, data->descC, msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase - single precision) .....");
    MORSE_strmm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, data->descC, data->descZ, msequence, mrequest);
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
    MORSE_MLE_szcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    MORSE_Sequence_Wait(data->sequence);
    VERBOSE(" Done.\n");
    free(streamdata);
    }
    */
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous - single precision)\n");
    VERBOSE("************************************************************\n");
}



double MORSE_smle_Tile(unsigned n, const double * theta, double * grad, void * MORSE_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- MORSE-sync-single precision.
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library. 
     * @param[in] MORSE_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0, zcpy_time=0.0;
    int N, NRHS, success;
    double flops =0.0;	

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
    START_TIMING(zcpy_time);
    if(data->iter_count==0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        MORSE_slacpy_Tile(MorseUpperLower, MORSE_descZ, MORSE_descZcpy); 
    else
    {	
        VERBOSE("Re-store the original Z vector...");
        MORSE_slacpy_Tile(MorseUpperLower, MORSE_descZcpy,MORSE_descZ);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(zcpy_time);	

    //****************************
    //   double *C = (double *) malloc(N * N * sizeof(double));

    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix (single precision)...");
    START_TIMING(matrix_gen_time);	
    //MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC, msequence, &mrequest[0], &data->l1, &data->l1, (double *)theta,  data->dm);    
    MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC, &data->l1, &data->l1, &data->lm, (double *)theta, data->dm, data->kernel_fun,  msequence, &mrequest[0]);
    STOP_TIMING(matrix_gen_time);
    MORSE_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");
    //double *C = (double *) malloc(N * N * sizeof(double));
    //       MORSE_Tile_to_Lapack( MORSE_descC, C, N);
    //   print_smatrix("test", 16, 16, C, 16);

    //double *C = (double *) malloc(N * N * sizeof(double));
    //MORSE_Tile_to_Lapack( MORSE_descC, C, N);
    //print_smatrix("testC", 16, 16, C, 16);

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma (single precision)...");
    START_TIMING(time_facto);
    success = MORSE_spotrf_Tile(MorseLower, MORSE_descC);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_SPOTRF(N);
    VERBOSE(" Done.\n");

    //       MORSE_Tile_to_Lapack( MORSE_descC, C, N);
    //       print_smatrix("test", 16, 16, C, 16);

    // MORSE_Tile_to_Lapack( MORSE_descC, C, N);
    //print_smatrix("testC", 16, 16, C, 16);
    //	exit(0);

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant (single precision) ...");
    START_TIMING(logdet_calculate);
    MORSE_MLE_smdet_Tile_Async(MORSE_descC, msequence, &mrequest[0], MORSE_descdet);
    MORSE_Sequence_Wait(msequence);
    //	printf("det: %f\n", data->det);
    logdet= 2 * (float)data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");
    //exit(0);

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system (single precision)...\n");
    START_TIMING(time_solve);
    MORSE_strsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_STRSM(MorseLeft,N, NRHS);
    VERBOSE(" Done.\n");    


    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function (single precision) ...");
    //dotp=0;
    //MORSE_MLE_core_ddotp_Async(MORSE_descZ,MORSE_descproduct,msequence, &mrequest[0]);
    //MORSE_Sequence_Wait(msequence);        

    MORSE_sgemm_Tile (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct); 

    loglik = -(N /2) + (N /2)*log (N) -(N / 2 ) * log(data->sdotp) -  0.5*logdet - (double) (N / 2.0) * log(2.0 * PI);
    VERBOSE(" Done.\n");

    data->variance = (1.0/N) * data->sdotp;

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(MORSE_My_Mpi_Rank() == 0)
    {
#endif
        //Print Iteration Summary
        //fprintf(stderr,"***************************************************\n");
        fprintf(stderr,"\n------ddotproduct: %.17g ", data->sdotp);

        fprintf(stderr,"\n------logdet: %.17g ", logdet);
        //      //fprintf(stderr,"------det: %.*e ", det);
        //      fprintf(stderr,"\n------expr2: %2.6f \n",((double) (N / 2) * log(2 * PI)));
        //fprintf(stderr," ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----LogLi: %2.6f\n", theta[0], theta[1], theta[2],loglik);
        //reformat
        printf(" %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %.17g\n", data->iter_count+1, data->variance, theta[1], theta[2],loglik);

        if(data->log == 1)
            fprintf(data->pFileLog, " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %.3g\n", data->iter_count+1,  data->variance, theta[1], theta[2],loglik);

        printf(" ---- Facto Time: %6.2f\n", time_facto);
        printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
        printf(" ---- dtrsm Time: %6.2f\n", time_solve);
        printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
        //fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", zcpy_time);
        printf(" ---- Total Time: %6.2f\n", matrix_gen_time+ time_facto + logdet_calculate + time_solve);
        //fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );
        printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
        //fprintf(stderr," ---- Peak Performance: %6.2f Gflops/s\n",  (ncores*p_grid*q_grid*16*2.3) );
        //fprintf(stderr,"***************************************************\n");
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}

double MORSE_smle_Tile_Async(unsigned n, const double * theta, double * grad, void * MORSE_data) {
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
    double loglik=0.0,  logdet=0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0, zcpy_time=0.0, flops = 0.0;
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
    START_TIMING(zcpy_time);
    if(data->iter_count == 0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        MORSE_slacpy_Tile_Async(MorseUpperLower, MORSE_descZ, MORSE_descZcpy, msequence, mrequest); 
    else
    {	
        VERBOSE("re-store the original Z vector (single precision)...");
        MORSE_slacpy_Tile_Async(MorseUpperLower, MORSE_descZcpy, MORSE_descZ, msequence, mrequest);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(zcpy_time);	


    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix (single precision)...");
    START_TIMING(matrix_gen_time);	
    //MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest, &data->l1, &data->l1, (double*) theta, data->dm);    
    MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC, &data->l1, &data->l1, &data->lm, (double *)theta, data->dm, data->kernel_fun,  msequence, &mrequest[0]);
    MORSE_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma (single precision)...");
    START_TIMING(time_facto);
    success = MORSE_spotrf_Tile_Async(MorseLower, MORSE_descC, msequence, mrequest);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_SPOTRF(N);
    VERBOSE(" Done.\n");

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant (single precision) ...");
    START_TIMING(logdet_calculate);
    MORSE_MLE_smdet_Tile_Async(MORSE_descC, msequence, &mrequest[0], MORSE_descdet);
    MORSE_Sequence_Wait(msequence);
    logdet= 2 * data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system (single precision) ...\n");
    START_TIMING(time_solve);
    MORSE_strsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ, msequence, mrequest);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_STRSM(MorseLeft, N, NRHS);
    VERBOSE(" Done.\n");    

    //Claculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function (single precision) ...");
    MORSE_sgemm_Tile_Async (MorseTrans, MorseNoTrans, 1, MORSE_descZ, MORSE_descZ, 0, MORSE_descproduct, msequence, mrequest); 
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

        fprintf(stderr," ---- Facto Time: %6.2f\n", time_facto);
        fprintf(stderr," ---- logdet Time: %6.2f\n", logdet_calculate);
        fprintf(stderr," ---- dtrsm Time: %6.2f\n", time_solve);
        fprintf(stderr," ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
        //fprintf(stderr," ---- re-store Z Vector Time: %6.2f\n", zcpy_time);
        fprintf(stderr," ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
        //fprintf(stderr," ---- Gflop (ignore): %6.2f\n", flops / 1e9 );    
        fprintf(stderr," ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
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


void MORSE_smle_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid, int mse_flag)
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
        return;
    }
    //Descriptors Creation
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZobs, NULL, MorseRealFloat, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
    if( mse_flag == 1)
    {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZactual, NULL, MorseRealFloat, dts, dts, dts * dts,  nZmiss, 1,  0, 0, nZmiss, 1, p_grid, q_grid);

        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descmse, &data->mserror, MorseRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZmiss, NULL, MorseRealFloat, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC12, NULL, MorseRealFloat, dts, dts, dts * dts, nZmiss, nZobs, 0, 0, nZmiss, nZobs, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC22, NULL, MorseRealFloat, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs, nZobs, p_grid, q_grid);

    //Initiate data descriptors
    data->descZmiss         = MORSE_descZmiss;
    data->descC12           = MORSE_descC12;
    data->descC22           = MORSE_descC22;
    data->descmse           = MORSE_descmse;
    data->descZactual       = MORSE_descZactual;
    data->descZobs          = MORSE_descZobs;
}


double MORSE_smle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n)
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
    //double *z = NULL, *streamdata = NULL;
    double time_solve = 0.0;
    double mat_gen_time = 0.0;
    double time_gemm = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;

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


    //Initiate data descriptors
    MORSE_descZmiss		= data->descZmiss;
    MORSE_descC12		= data->descC12;
    MORSE_descC22		= data->descC22;
    MORSE_descmse		= data->descmse;
    MORSE_descZactual	= data->descZactual;
    MORSE_descZobs		= data->descZobs;

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
    theta[0]=data->variance;

    printf("estimated parameters: %f - %f - %f\n", theta[0], theta[1], theta[2]);
    //Copy data to vectors 
    VERBOSE("Copy measurments vector to descZobs descriptor...");
    //MORSE_MLE_szcpy_Tile_Async(MORSE_descZobs, Zobs, msequence, mrequest);
    MORSE_Lapack_to_Tile( Zobs, nZobs, MORSE_descZobs);
    VERBOSE(" Done.\n");

    if( Zactual != NULL)	
    {
        VERBOSE("Copy actual measurments vector to descZactual descriptor...");
        //MORSE_MLE_szcpy_Tile_Async(MORSE_descZactual, Zactual, msequence, mrequest);
        MORSE_Lapack_to_Tile( Zactual, nZmiss, MORSE_descZactual);
        VERBOSE(" Done.\n");
    }

    MORSE_Sequence_Wait(msequence);



    START_TIMING(mat_gen_time);
    //Generate C22 covariance matrix
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage - single precision)");
    MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZobs);
    VERBOSE(" Done.\n");


    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage - single precision)");
    //MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  &data->lmiss, &data->lobs, theta, data->dm);
    MORSE_MLE_dcmg_Tile_Async(MorseLower, MORSE_descC12, &data->lmiss, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZmiss);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);



    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage - single precision)");
    MORSE_sposv_Tile(MorseLower, MORSE_descC22, MORSE_descZobs);
    flops = flops + FLOPS_SPOTRF(nZobs);
    flops = flops + FLOPS_STRSM(MorseLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);


    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage - single precision)");
    MORSE_sgemm_Tile (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss);
    flops = flops + FLOPS_SGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);


    //return back descZmiss to zmiss vector
    MORSE_Tile_to_Lapack( MORSE_descZmiss, Zmiss, nZmiss);

    //Estimate Mean Square Error
    if( Zactual != NULL)
    {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
        MORSE_MLE_smse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
        MORSE_Sequence_Wait(msequence);
        VERBOSE(" Done.\n");	
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

        write_prediction_result("predict_result.dat", n, data->hicma_acc, 0, 0, data->mserror, (mat_gen_time+time_solve+ time_gemm), (flops / 1e9 / (time_solve)));

#if defined(CHAMELEON_USE_MPI)
    }
#endif


    return data->mserror;

}



double MORSE_smle_Predict_Tile_Async(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n)
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
    double time_solve = 0.0;
    double mat_gen_time = 0.0;
    double time_gemm = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;

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

    //Initiate data descriptors.
    MORSE_descZmiss         = data->descZmiss;
    MORSE_descC12           = data->descC12;
    MORSE_descC22           = data->descC22;
    MORSE_descmse           = data->descmse;
    MORSE_descZactual       = data->descZactual;
    MORSE_descZobs          = data->descZobs;

    //Copy data to vectors.
    VERBOSE("Copy measurments vector to descZobs descriptor (single precision)...");
    //MORSE_Lapack_to_Tile_Async( Zobs, nZobs, MORSE_descZobs, msequence, mrequest);
    VERBOSE(" Done.\n");

    if( Zactual != NULL)
    {
        VERBOSE("Copy actual measurments vector to descZactual descriptor (single precision)...");
        //      MORSE_Lapack_to_Tile_Async( Zactual, nZmiss, MORSE_descZactual, msequence, mrequest);
        VERBOSE(" Done.\n");
    }


    START_TIMING(mat_gen_time);
    //Generate C22 covariance matrix.
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage - single precision)");
    //MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC22, msequence, mrequest,  &data->lobs, &data->lobs, theta, data->dm);
    MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
    VERBOSE(" Done.\n");


    //Generate C12 covariance matrix.
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage - single precision)");
    //MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC12, msequence, mrequest,  &data->lmiss, &data->lobs, theta, data->dm);
    MORSE_MLE_scmg_Tile_Async(MorseLower, MORSE_descC12, &data->lmiss, &data->lobs, &data->lm, theta, data->dm, data->kernel_fun,  msequence, mrequest);
    MORSE_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);



    START_TIMING(time_solve);
    //Start prediction.
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage - single precision)");
    MORSE_sposv_Tile_Async(MorseLower, MORSE_descC22, MORSE_descZobs, msequence, mrequest);
    flops = flops + FLOPS_SPOTRF(nZobs);
    flops = flops + FLOPS_STRSM(MorseLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);


    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage - single precision)");
    MORSE_sgemm_Tile_Async (MorseNoTrans, MorseNoTrans, 1, MORSE_descC12, MORSE_descZobs, 0, MORSE_descZmiss, msequence, mrequest);
    flops = flops + FLOPS_SGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);


    //return back descZmiss to zmiss vector
    //MORSE_Tile_to_Lapack_Async( MORSE_descZmiss, Zmiss, nZmiss, msequence, mrequest);

    //Estimate Mean Square Error
    if( Zactual != NULL)
    {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage - single precision) \n");
        MORSE_MLE_smse_Tile_Async(MORSE_descZactual, MORSE_descZmiss, MORSE_descmse, msequence, mrequest);
        VERBOSE(" Done.\n");
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

        write_prediction_result("predict_result.dat", n, data->hicma_acc, 0, 0, data->mserror, (mat_gen_time+time_solve+ time_gemm), (flops / 1e9 / (time_solve)));
#if defined(CHAMELEON_USE_MPI)
    }
#endif


    return data->mserror;

}



//init Chameleon descriptors
void MORSE_smle_Call(MLE_data  *data, int ncores,int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
    //! //Initiate MORSE and allocate different descriptors for
    /*!  CHAMELEON
     * Returns MLE_data data with initial values and new descriptors locations.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] ncores: number of CPU workers.
     * @param[in] gpus: number of GPU workers.
     * @param[in] dts: tile size (MB).
     * @param[in] p_grid: p_grid.
     * @param[in] q_grid: q_grid.
     * @param[in] N: number of spatial locations.	
     * @param[in] nZmiss: number of missing values (unknown observations).
     * @param[in] nZobs: number of observed values (known observations).
     * */
{

    MORSE_sequence_t *msequence;
    MORSE_request_t mrequest[2] = { MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER };
    MORSE_desc_t *MORSE_descC	= NULL;
    MORSE_desc_t *MORSE_descZ	= NULL;
    MORSE_desc_t *MORSE_descZcpy	= NULL;
    MORSE_desc_t *MORSE_descproduct	= NULL;
    MORSE_desc_t *MORSE_descdet	= NULL;

    // For ditributed system and should be removed
    float *Zcpy = (float *) malloc(N * sizeof(float));

    //Identifies a set of routines sharing common exception handling.
    MORSE_Sequence_Create(&msequence);

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descC, NULL , MorseRealFloat, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZ, NULL, MorseRealFloat, dts, dts, dts * dts, N, 1,  0, 0, N, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descZcpy, Zcpy, MorseRealFloat, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descproduct, &data->sdotp, MorseRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&MORSE_descdet, &data->det, MorseRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    //Fill data struct
    data->descC		= MORSE_descC;
    data->descZ		= MORSE_descZ;
    data->descZcpy		= MORSE_descZcpy;
    data->descdet		= MORSE_descdet;
    data->descproduct	= MORSE_descproduct;
    data->sequence		= msequence;
    data->request		= mrequest;
    //stop gsl error handler
    gsl_set_error_handler_off () ;

}


