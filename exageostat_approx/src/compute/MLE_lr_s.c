/**
 *
 * @file MLE.c
 *
 *  
 *  ExaGeoStat is a software package provided by KAUST,
 *  King Abdullah University of Science and Technology - ECRC
 *
 * @version 1.2.0
 * @author Sameh Abdulah
 * @date 2022-11-09
 * @generated d Fri Nov 22 15:11:13 2016
 *
 **/
#include "../include/MLE_lr_s.h"

//***************************************************************************************
//HiCMA global variables.
extern int store_only_diagonal_tiles = 1;
extern int print_index = 0;
extern STARSH_blrf *mpiF;
extern int use_scratch = 1;
extern int global_check = 0;  //used to create dense matrix for accuracy check
extern int print_mat = 0;

//Generate Observations Vector (Z) for testing Maximum Likelihood function -- CHAM-sync
int HICMA_MLE_szvg_Tile(MLE_data *data, double* Nrand, double* initial_theta, int n, int dts, int log, int p_grid,
                        int q_grid) {

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
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC, msequence, mrequest, &data->l1, &data->l1, initial_theta,
                                   data->dm);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("LR: Generate Normal Random Distribution Vector Z (Synthetic Dataset Dense Generation Phase) .....");
    EXAGEOSTAT_MLE_szcpy_Tile_Async(CHAMELEON_descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("LR: Cholesky factorization of Sigma (Synthetic Dataset Dense Generation Phase) .....");
    int success = CHAMELEON_spotrf_Tile(ChamLower, CHAM_descC);
    SUCCESS(success, "LR: Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("LR: Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Dense Generation Phase) .....");
    CHAMELEON_strmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAMELEON_descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if (log == 1) {
        double* z;
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

void EXAGEOSTAT_TLR_MLE_zcpy(MLE_data *data, double* streamdata) {
    RUNTIME_sequence_t *hsequence = (RUNTIME_sequence_t *) data->hsequence;
    RUNTIME_request_t *hrequest = (RUNTIME_request_t *) data->hrequest;
    VERBOSE("LR: Copy Z from vector to decriptor.\n");
    EXAGEOSTAT_MLE_dzcpy_Tile_Async(data->hicma_descZ, streamdata, hsequence, hrequest);
    CHAMELEON_Sequence_Wait(hsequence);
    VERBOSE("LR: Done Z copying step.\n");
    VERBOSE("************************************************************\n");
}

int HICMA_MLE_szvg_Tile_Async(MLE_data *data, double* Nrand, double* initial_theta, int n, int dts, int log, int p_grid,
                              int q_grid) {
    //// TODO: Missing Implementation.
}

//compute MLE function
double HICMA_smle_Tile(unsigned n, const double* theta, double* grad, void *app_data) {

    //Initialization
    double loglik = 0.0;
    int compress_diag = 0;
    //double  ddotproduct	= 0.0,
    double logdet = 0.0;
    double time_facto = 0.0,
            time_solve = 0.0,
            logdet_calculate = 0.0,
            matrix_gen_time = 0.0,
            test_time = 0.0;
    double flops = 0.0;
    int success, maxrank, acc;
    //int verbose;
    int N, NRHS;
    int lts;//, dts;
    int hicma_data_type;
    int i = 0;

    MLE_data *data = (MLE_data *) app_data;
    CHAM_desc_t *hicma_descC = (CHAM_desc_t *) data->hicma_descC;
    CHAM_desc_t *hicma_descCD = (CHAM_desc_t *) data->hicma_descCD;
    CHAM_desc_t *hicma_descCUV = (CHAM_desc_t *) data->hicma_descCUV;
    CHAM_desc_t *hicma_descCrk = (CHAM_desc_t *) data->hicma_descCrk;
    CHAM_desc_t *hicma_descZ = (CHAM_desc_t *) data->hicma_descZ;
    CHAM_desc_t *hicma_descZcpy = (CHAM_desc_t *) data->hicma_descZcpy;
    CHAM_desc_t *hicma_descdet = (CHAM_desc_t *) data->hicma_descdet;
    CHAM_desc_t *hicma_descproduct = (CHAM_desc_t *) data->hicma_descproduct;
    RUNTIME_sequence_t *hsequence = (RUNTIME_sequence_t *) data->hsequence;
    RUNTIME_request_t *hrequest = (RUNTIME_request_t *) data->hrequest;
    //verbose 	  		= data->verbose;
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
            double* z = (double* ) malloc(N * sizeof(double));
            CHAMELEON_Tile_to_Lapack(hicma_descZ, z, N);
            CHAMELEON_Lapack_to_Tile(z, N, hicma_descZ);
            free(z);
            CHAMELEON_Desc_Destroy(&hicma_descZ);
        }
        CHAMELEON_slacpy_Tile(ChamUpperLower, hicma_descZ, hicma_descZcpy);
    }

    //Matrix generation part.       
    VERBOSE("LR:Generate New Covariance Matrix (single precision)...");
    START_TIMING(matrix_gen_time);

    HICMA_problem_t hicma_problem;
    hicma_problem.theta = (double* ) theta;
    hicma_problem.noise = 1e-4;
    hicma_problem.ndim = 2;


    // I need to change the location struct to vector array (TO DO)
    double* xycoord = (double* ) malloc(2 * N * sizeof(double));
    for (i = 0; i < N; i++) {
        xycoord[i] = data->l1.x[i];
        xycoord[N + i] = data->l1.y[i];
    }


    if (strcmp(data->dm, "gc") == 0)
        printf("gcd metric is used: %d\n", STARSH_SPATIAL_MATERN2_GCD);
    else
        printf("ed metric is used: %d\n", STARSH_SPATIAL_MATERN2_SIMD);

    hicma_problem.point = xycoord;
    hicma_problem.kernel_type = strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_MATERN2_GCD : STARSH_SPATIAL_MATERN2_SIMD;
    HICMA_zgenerate_problem(hicma_data_type, 'S', 0, N, lts, hicma_descCUV->mt, hicma_descCUV->nt, &hicma_problem);
    mpiF = hicma_problem.starsh_format;
    HICMA_zgytlr_Tile(ChamLower, hicma_descCUV, hicma_descCD, hicma_descCrk, 0, maxrank, pow(10, -1.0 * acc),
                      compress_diag, hicma_descC);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");

    //******************************
    VERBOSE("LR: re-Copy z (single precision)...");
    START_TIMING(test_time);
    //re-store old Z
    CHAMELEON_slacpy_Tile(ChamUpperLower, hicma_descZcpy, hicma_descZ);
    STOP_TIMING(test_time);
    VERBOSE(" Done.\n");

    //**************
    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("LR: Cholesky factorization of Sigma (single precision)...");
    START_TIMING(time_facto);

    success = HICMA_zpotrf_Tile(ChamLower, hicma_descCUV, hicma_descCD, hicma_descCrk, 0, maxrank, pow(10, -1.0 * acc));

    SUCCESS(success, "LR: Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    STOP_TIMING(time_facto);
    flops = flops + FLOPS_DPOTRF(N);
    VERBOSE(" Done.");

    VERBOSE("LR:Calculating the log determinant (single precision)...");
    START_TIMING(logdet_calculate);
    data->det = 0;
    //// TODO: RESOLVE THIS!
    HICMA_MLE_smdet_Tile_Async(hicma_descCD, hsequence, &hrequest[0], hicma_descdet);
    CHAMELEON_Sequence_Wait(hsequence);
    logdet = 2 * data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("LR:Solving the linear system (single precision)...");
    START_TIMING(time_solve);
    //Compute triangular solve LC*X = Z

    HICMA_ztrsmd_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, hicma_descCUV, hicma_descCD, hicma_descCrk,
                      hicma_descZ, maxrank);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_DTRSM(ChamLeft, N, NRHS);
    VERBOSE(" Done.");

    VERBOSE("LR:Calculating dot product (single precision)...");
    CHAMELEON_sgemm_Tile(ChamTrans, ChamNoTrans, 1, hicma_descZ, hicma_descZ, 0, hicma_descproduct);
    loglik = -0.5 * ((MLE_data *) data)->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);

    VERBOSE(" Done.");

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    //Print Iteration Summary

    printf(" %3d- Model Parameters (varinace, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n",
           ((MLE_data *) data)->iter_count + 1, theta[0], theta[1], theta[2], loglik);

    if (((MLE_data *) data)->log == 1)
        fprintf(((MLE_data *) data)->pFileLog,
                " %3d- Model Parameters (varinace, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n",
                ((MLE_data *) data)->iter_count + 1, theta[0], theta[1], theta[2], loglik);
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    ((MLE_data *) data)->iter_count++;
    // for experiments
    ((MLE_data *) data)->avg_exec_time_per_iter += matrix_gen_time + time_facto + logdet_calculate + time_solve;
    ((MLE_data *) data)->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    ((MLE_data *) data)->final_loglik = loglik;
    return loglik;
}

double HICMA_smle_Tile_Async(unsigned n, const double* theta, double* grad, void *data) {

    //Initialization
}

void
HICMA_smle_Predict_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs, int lts, int p_grid, int q_grid, int mse_flag) {
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
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZobs, NULL, ChamRealFloat, lts, lts, lts * lts, nZobs, 1, 0, 0,
                                    nZobs, 1, p_grid, q_grid);

    if (mse_flag == 1) {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZactual, NULL, ChamRealFloat, lts, lts, lts * lts, nZmiss, 1, 0,
                                        0, nZmiss, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealFloat, lts, lts, lts * lts, 1, 1, 0, 0,
                                        1, 1, p_grid, q_grid);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZmiss, NULL, ChamRealFloat, lts, lts, lts * lts, nZmiss, 1, 0, 0,
                                    nZmiss, 1, p_grid, q_grid);


    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC12, NULL, ChamRealFloat, lts, lts, lts * lts, nZmiss, nZobs, 0, 0,
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

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC, NULL, ChamRealFloat, MBC, NBC, MBC * NBC, MC, NC, 0, 0, MC, NC,
                                    p_grid, q_grid);

    MBD = lts;
    NBD = lts;
    MD = nZobs;
    ND = MBD;

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22D, NULL, ChamRealFloat, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND,
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
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22UV, NULL, ChamRealFloat, MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0, MUV,
                                    NUV, p_grid, q_grid);

    //CUV Descriptor
    MBrk = 1;
    NBrk = 1;
    Mrk = hicma_descC22UV->mt;
    Nrk = hicma_descC22UV->mt;
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22rk, NULL, ChamRealFloat, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk,
                                    Nrk, p_grid, q_grid);

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
double HICMA_smle_Predict_Tile(MLE_data *CHAM_data, double* theta, int nZmiss, int nZobs, double* Zobs, double* Zactual,
                               double* Zmiss, int n, int lts)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-sync
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

    VERBOSE(" LR: Generate C22 Covariance Matrix... (Prediction Stage)");
    START_TIMING(mat_gen_time);

    HICMA_problem_t hicma_problem;
    hicma_problem.theta = (double* ) theta;
    hicma_problem.noise = 1e-4;
    hicma_problem.ndim = 2;

    // I need to change the location struct to vector array (TO DO)
    double* xycoord = (double* ) malloc(2 * nZobs * sizeof(double));
    int i;
    for (i = 0; i < nZobs; i++) {
        xycoord[i] = data->lobs.x[i];
        xycoord[nZobs + i] = data->lobs.y[i];
    }

    hicma_problem.point = xycoord;
    hicma_problem.kernel_type = strcmp(data->dm, "gc") == 0 ? STARSH_SPATIAL_MATERN2_GCD : STARSH_SPATIAL_MATERN2_SIMD;
    HICMA_zgenerate_problem(hicma_data_type, 'S', 0, nZobs, lts, hicma_descC22UV->mt, hicma_descC22UV->nt,
                            &hicma_problem);
    mpiF = hicma_problem.starsh_format;
    HICMA_zgytlr_Tile(ChamLower, hicma_descC22UV, hicma_descC22D, hicma_descC22rk, 0, maxrank, pow(10, -1.0 * acc),
                      compress_diag, hicma_descC);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix
    VERBOSE("LR: Generate C12 Covariance Matrix... (Prediction Stage)");
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC12, hsequence, hrequest, &data->lmiss, &data->lobs, theta,
                                   data->dm);
    CHAMELEON_Sequence_Wait(hsequence);
    //flops = flops + FLOPS_DPOTRF(nZmiss);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);

    //***************************************
    START_TIMING(time_solve);
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
    HICMA_zpotrf_Tile(ChamLower, hicma_descC22UV, hicma_descC22D, hicma_descC22rk, 0, maxrank, pow(10, -1.0 * acc));
    HICMA_ztrsmd_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, hicma_descC22UV, hicma_descC22D,
                      hicma_descC22rk, CHAMELEON_descZobs, maxrank);
    HICMA_ztrsmd_Tile(ChamLeft, ChamLower, ChamTrans, ChamNonUnit, 1, hicma_descC22UV, hicma_descC22D, hicma_descC22rk,
                      CHAMELEON_descZobs, maxrank);
    flops = flops + FLOPS_DPOTRF(nZobs);
    flops = flops + FLOPS_DTRSM(ChamLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);

    //**********************************dgemm
    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");
    CHAMELEON_sgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAMELEON_descZobs, 0, CHAMELEON_descZmiss);
    flops = flops + FLOPS_DGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //return back descZmiss to zmiss vector
    CHAMELEON_Tile_to_Lapack(CHAMELEON_descZmiss, Zmiss, nZmiss);

    //Estimate Mean Square Error
    if (Zactual != NULL) {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
        EXAGEOSTAT_MLE_smse_Tile_Async(CHAMELEON_descZactual, CHAMELEON_descZmiss, CHAM_descmse, hsequence, hrequest);
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
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n",
                nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

    write_prediction_result("predict_result.dat", n, data->hicma_acc, data->mserror,
                            (mat_gen_time + time_solve + time_gemm), (flops / 1e9 / (time_solve)));

#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;
}

//init Hicma decriptors
void
HICMA_smle_Call(MLE_data *data, int ncores, int gpus, int lts, int p_grid, int q_grid, int N, int nZobs, int nZmiss) {

    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAM_desc_t *hicma_descC = NULL;
    CHAM_desc_t *hicma_descZ = NULL;
    CHAM_desc_t *hicma_descZcpy = NULL;
    CHAM_desc_t *hicma_descproduct = NULL;
    CHAM_desc_t *hicma_descdet = NULL;
    CHAM_desc_t *CHAMELEON_descZmiss = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAMELEON_descZactual = NULL;
    CHAM_desc_t *CHAMELEON_descZobs = NULL;
    CHAM_desc_t *hicma_descCD = NULL;
    CHAM_desc_t *hicma_descCUV = NULL;
    CHAM_desc_t *hicma_descCrk = NULL;
    CHAM_desc_t *hicma_descC12D = NULL;
    CHAM_desc_t *hicma_descC12UV = NULL;
    CHAM_desc_t *hicma_descC12rk = NULL;
    CHAM_desc_t *hicma_descC22D = NULL;
    CHAM_desc_t *hicma_descC22UV = NULL;
    CHAM_desc_t *hicma_descC22rk = NULL;
    int MBC, NBC, MC, NC;
    int MBD, NBD, MD, ND;
    int MBUV, NBUV, MUV, NUV;
    int MBrk, NBrk, Mrk, Nrk;

    //For ditributed system and should be removed
    double* Zcpy = (double* ) malloc(N * sizeof(double));

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

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC, NULL, ChamRealFloat, MBC, NBC, MBC * NBC, MC, NC, 0, 0, MC, NC,
                                    p_grid, q_grid);
    printf("(1)%d - %d\n", MC, NC);

    //CAD Descriptor
    MBD = lts;
    NBD = lts;
    MD = N;
    ND = MBD;
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descCD, NULL, ChamRealFloat, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND,
                                    p_grid, q_grid);
    printf("(2)%d - %d\n", MD, ND);
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
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descCUV, NULL, ChamRealFloat, MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0, MUV,
                                    NUV, p_grid, q_grid);
    printf("(3)%d - %d\n", MUV, NUV);

    //CUV Descriptor
    MBrk = 1;
    NBrk = 1;
    Mrk = hicma_descCUV->mt;
    Nrk = hicma_descCUV->mt;
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descCrk, NULL, ChamRealFloat, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk,
                                    Nrk, p_grid, q_grid);
    CHAMELEON_Sequence_Create(&msequence);
    printf("(4)%d - %d\n", Mrk, Nrk);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descZ, NULL, ChamRealFloat, lts, lts, lts * lts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descZcpy, Zcpy, ChamRealFloat, lts, lts, lts * lts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descproduct, &data->dotp, ChamRealFloat, lts, lts, lts * lts, 1, 1, 0, 0, 1,
                                    1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descdet, &data->det, ChamRealFloat, lts, lts, lts * lts, 1, 1, 0, 0, 1, 1,
                                    p_grid, q_grid);


    if (nZmiss != 0) {
        if (strcmp(data->actualZFPath, "") == 0) {

            EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZobs, &Zcpy[nZmiss], ChamRealFloat, lts, lts, lts * lts,
                                            nZobs, 1, 0, 0, nZobs, 1, p_grid, q_grid);
            EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZactual, Zcpy, ChamRealFloat, lts, lts, lts * lts, nZmiss, 1,
                                            0, 0, nZmiss, 1, p_grid, q_grid);
        } else {
            CHAMELEON_descZobs = hicma_descZcpy;
            EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZactual, NULL, ChamRealFloat, lts, lts, lts * lts, nZmiss, 1,
                                            0, 0, nZmiss, 1, p_grid, q_grid);
        }


        //C12AD Descriptor    
        MBD = lts;
        NBD = lts;
        MD = nZmiss;
        ND = MBD;
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC12D, NULL, ChamRealFloat, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND,
                                        p_grid, q_grid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * data->hicma_maxrank;
        MUV = nZmiss;
        NUV = 2 * MUV / lts * data->hicma_maxrank;
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC12UV, NULL, ChamRealFloat, MBUV, NBUV, MBUV * NBUV, MBUV, NBUV, 0,
                                        0, MBUV, NBUV, p_grid, q_grid);

        //C12Ark Descriptor
        MBrk = 1;
        NBrk = 1;
        Mrk = hicma_descC12UV->mt;
        Nrk = hicma_descC12UV->mt;
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC12rk, NULL, ChamRealFloat, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,
                                        Mrk, Nrk, p_grid, q_grid);

        //C11AD Descriptor
        MBD = lts;
        NBD = lts;
        MD = nZobs;
        ND = MBD;
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22D, NULL, ChamRealFloat, MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND,
                                        p_grid, q_grid);

        //C12UV Descriptor
        MBUV = lts;
        NBUV = 2 * data->hicma_maxrank;
        MUV = nZobs;
        NUV = 2 * MUV / lts * data->hicma_maxrank;
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22UV, NULL, ChamRealFloat, MBUV, NBUV, MBUV * NBUV, MBUV, NBUV, 0,
                                        0, MBUV, NBUV, p_grid, q_grid);
        //C12Ark Descriptor            
        MBrk = 1;
        NBrk = 1;
        Mrk = hicma_descC22UV->mt;
        Nrk = hicma_descC22UV->mt;
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&hicma_descC22rk, NULL, ChamRealFloat, MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0,
                                        Mrk, Nrk, p_grid, q_grid);

        //Other descriptors
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAMELEON_descZmiss, NULL, ChamRealFloat, lts, lts, lts * lts, nZmiss, 1, 0, 0,
                                        nZmiss, 1, p_grid, q_grid);
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealFloat, lts, lts, lts * lts, 1, 1, 0, 0,
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
    data->hicma_descproduct = hicma_descproduct;
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
    data->hsequence = msequence;
    data->hrequest = mrequest;
    data->mserror = 0;
    //stop gsl error handler
    gsl_set_error_handler_off();
}