/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file rwrappers.c
 *
 * ExaGeoStat R-wrapper functions.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/rwrappers.h"
#include "MLE_lr.h"
#include "hicma/include/hicma.h"
#include "hicma/hicma_ext/control/hicma_context.h"


static int get_num_params(char *kernel_fun) {
    int num_params = 0;

    if (strcmp(kernel_fun, "univariate_matern_stationary") == 0)
        num_params = 3;
    else if (strcmp(kernel_fun, "univariate_matern_non_stationary") == 0)
        num_params = 9;
    else if (strcmp(kernel_fun, "bivariate_matern_flexible") == 0)
        num_params = 13;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious") == 0)
        num_params = 6;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious2") == 0)
        num_params = 6;
    else if (strcmp(kernel_fun, "trivariate_matern_parsimonious") == 0)
        num_params = 10;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        num_params = 6;
    else if (strcmp(kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        num_params = 7;
    else if (strcmp(kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        num_params = 10;
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    return num_params;
}

void gen_z_givenlocs_exact(double* x, int *xlen, double* y,
                           int *ylen, int *kernel, double* theta,
                           int *thetalen, int *dmetric, int *n,
                           int *ncores, int *gpus, int *ts,
                           int *p_grid, int *q_grid, int *veclen,
                           double* globalvec)
//! direct function to generate synthetics datasets (X, Y) 2D locations Z measurement vector.
/*!  -- using dense or approximate computation
 * Returns Z observation vector
 * @param[in] x:                Pointer to the x-dim vector.
 * @param[in] xlen:             Pointer to the length of x-dim vector.
 * @param[in] y:                Pointer to the y-dim vector.
 * @param[in] ylen:             Pointer to the length of y-dim vector.
 * @param[in] kernel:           Pointer to the computation kernel.
 * @param[in] theta:            Pointer to the theta vector.
 * @param[in] thetalen:         Pointer to the length of y-dim vector.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] veclen:           Pointer to the length of global vec (R memory space).
 * @param[in] globalvec:        Pointer to R memory space (x:y:z).
 * */
{
    //Initialization
    int j = 0, log = 0, verbose = 0, async = 0, comp_mode = 0, metric_mode = 0;
    double* initial_theta, *localvec;
    MLE_data data;
    location *locations;
    int iseed[4] = {0, 0, 0, 1};
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAM_descZ = NULL;
    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    //initialize globalvec output
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = 0;

    //Memory allocation
    locations = (location *) malloc(sizeof(location));
    initial_theta = (double* ) malloc((*thetalen) * sizeof(double));
    localvec = (double* ) malloc((*veclen) * sizeof(double));

    //set XY locs
    locations->x = x;
    locations->y = y;

    //Set data struct values based on inputs
    comp_mode = 0;  //Usually exact
    data.l1 = *locations;
    if (*kernel == 0)
        data.kernel_fun = "univariate_matern_stationary";
    else if (*kernel == 1)
        data.kernel_fun = "univariate_matern_nuggets_stationary";
    else if (*kernel == 2)
        data.kernel_fun = "bivariate_matern_flexible";
    else if (*kernel == 3)
        data.kernel_fun = "bivariate_matern_parsimonious";
    else if (*kernel == 4)
        data.kernel_fun = "trivariate_matern_parsimonious";
    else if (*kernel == 5)
        data.kernel_fun = "univariate_spacetime_matern_stationary";
    else if (*kernel == 6)
        data.kernel_fun = "bivariate_spacetime_matern_stationary";
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    data.computation = comp_mode == 0 ? "exact" : "appro";
    data.dm = metric_mode == 0 ? "ed" : "gcd";
    data.verbose = verbose;
    data.async = async;
    data.obsFPath = "";
    data.log = log;
    gsl_set_error_handler_off();

    //Assign initial_ theta vector to generate synthetic dataset
    for (j = 0; j < *thetalen; j++)
        initial_theta[j] = theta[j];

    //Nomral random generation of e -- ei~N(0, 1) to generate Z
    LAPACKE_dlarnv(3, iseed, *n, &localvec[0]);

    //Create Descriptors
    CHAMELEON_Sequence_Create(&msequence);
    CHAMELEON_Desc_Create(&CHAM_descC, NULL, ChamRealDouble, *ts, *ts, *ts * *ts, *n, *n, 0, 0, *n, *n, *p_grid,
                          *q_grid);
    CHAMELEON_Desc_Create(&CHAM_descZ, NULL, ChamRealDouble, *ts, *ts, *ts * *ts, *n, 1, 0, 0, *n, 1, *p_grid, *q_grid);
    data.sequence = msequence;
    data.request = mrequest;
    data.descC = CHAM_descC;
    data.descZ = CHAM_descZ;

    //Main algorithm call
    MLE_zvg(&data, localvec, initial_theta, *n, *ts, log, *p_grid, *q_grid);

    CHAMELEON_Tile_to_Lapack(data.descZ, localvec, *n);
    //copy local vector to global vector in R memory space
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = *(localvec + j);

    //Destory descriptors
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descC);
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descZ);

    //free memory
    free(initial_theta);
    free(localvec);
}

void gen_z_exact(int *kernel, double* theta, int *thetalen,
                 int *dmetric, int *n, int *seed,
                 int *ncores, int *gpus, int *ts,
                 int *p_grid, int *q_grid, int *veclen,
                 double* globalvec)
//! direct function to generate synthetics datasets (X, Y) 2D locations Z measurement vector.
/*!  -- using dense or approximate computation
 * Returns Z observation vector
 * @param[in] kernel:           Pointer to the used kernel
 * @param[in] theta:            Pointer to the theta vector.
 * @param[in] thetalen:         Pointer to the length of y-dim vector.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] seed:             Pointer to the random seed.
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] veclen:           Pointer to the length of global vec (R memory space).
 * @param[in] globalvec:        Pointer to R memory space (x:y:z).
 * */
{

    //Initialization
    int i = 0, j = 0, log = 0, verbose = 0, async = 0, comp_mode = 0, metric_mode = 0;
    double* initial_theta, *localvec;
    MLE_data data;
    location *locations;
    int iseed[4] = {*seed, *seed, *seed, 1};
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAM_descZ = NULL;
    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    //initialize globalvec
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = 0;

    //Memory allocation
    initial_theta = (double* ) malloc(*thetalen * sizeof(double));
    localvec = (double* ) malloc(*veclen * sizeof(double));

    //Generate XY locs
    locations = GenerateXYLoc(*n, *seed);
    data.l1 = *locations;
    //Assign x, y vectors to local space (C memeory space)
    for (i = 0; i < *n; i++)
        localvec[i] = data.l1.x[i];
    for (i = *n; i < 2 * (*n); i++)
        localvec[i] = data.l1.y[i - (*n)];
    //Set data struct values based on inputs
    comp_mode = 0;  //Usually exact

    if (*kernel == 0)
        data.kernel_fun = "univariate_matern_stationary";
    else if (*kernel == 1)
        data.kernel_fun = "univariate_matern_nuggets_stationary";
    else if (*kernel == 2)
        data.kernel_fun = "bivariate_matern_flexible";
    else if (*kernel == 3)
        data.kernel_fun = "bivariate_matern_parsimonious";
    else if (*kernel == 4)
        data.kernel_fun = "trivariate_matern_parsimonious";
    else if (*kernel == 5)
        data.kernel_fun = "univariate_spacetime_matern_stationary";
    else if (*kernel == 6)
        data.kernel_fun = "bivariate_spacetime_matern_stationary";
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    data.l1 = *locations;
    data.computation = comp_mode == 0 ? "exact" : "appro";
    data.dm = metric_mode == 0 ? "ed" : "gcd";
    data.verbose = verbose;
    data.async = async;
    data.obsFPath = "";
    data.log = log;
    data.precision = 0;
    gsl_set_error_handler_off();

    //Assign initial_ theta vector to generate synthetic dataset
    for (j = 0; j < *thetalen; j++)
        initial_theta[j] = theta[j];

    //Nomral random generation of e -- ei~N(0, 1) to generate Z
    LAPACKE_dlarnv(3, iseed, *n, &localvec[2 * *n]);

    //Create Descriptors
    CHAMELEON_Sequence_Create(&msequence);
    CHAMELEON_Desc_Create(&CHAM_descC, NULL, ChamRealDouble, *ts, *ts, *ts * *ts, *n, *n, 0, 0, *n, *n, *p_grid,
                          *q_grid);
    CHAMELEON_Desc_Create(&CHAM_descZ, NULL, ChamRealDouble, *ts, *ts, *ts * *ts, *n, 1, 0, 0, *n, 1, *p_grid, *q_grid);
    data.sequence = msequence;
    data.request = mrequest;
    data.descC = CHAM_descC;
    data.descZ = CHAM_descZ;

    //Main algorithm call
    MLE_zvg(&data, &localvec[2 * *n], initial_theta, *n, *ts, log, *p_grid, *q_grid);

    CHAMELEON_Tile_to_Lapack(data.descZ, &localvec[2 * *n], *n);
    //copy local vector to global vector in R memory space
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = *(localvec + j);

    //Destory descriptors
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descC);
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descZ);

    //free memory
    free(initial_theta);
    free(localvec);
}

static void mle_general_tlr(int *kernel, double* x, int *xlen, double* y,
                            int *ylen, double* z, int *zlen,
                            double* clb, int *clblen, double* cub,
                            int *cublen, int *computation, int *diag_thick,
                            int *lr_acc, int *lr_maxrank, int *dmetric,
                            int *n, double* opt_tol, int *opt_max_iters,
                            int *ncores, int *gpus, int *dts, int *lts,
                            int *p_grid, int *q_grid,
                            double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] kernel_fun:       Pointer to the stationary/nono-stationary kernel.
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector).
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] diag_thick:       Pointer to the used diagonal thick in the case of DST computation..
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] opt_tol:        Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
 * */
{
    //initialization
    int i = 0, j = 0;
    //double time_opt = 0.0, max_loglik = 0.0;
    //Log is 0 in R
    int log = 0;
    //Verbose is 0 in R
    int verbose = 0;
    //Async is 0 in R
    int async = 0;
    int comp_mode = *computation;
    int metric_mode = *dmetric;
    double opt_f;
    nlopt_opt opt;
    MLE_data data;
    double* starting_theta;

    HICMA_context_t *hicmatxt;
    hicmatxt = hicma_context_self();
    if(hicmatxt == NULL){
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }
    printf("Done\n");

    //Initialize the data struct
    init_data_values(&data);
    int num_params = 3;
    if (*kernel == 0) {
        data.kernel_fun = "univariate_matern_stationary";
        num_params = 3;
    } else if (*kernel == 1) {
        data.kernel_fun = "univariate_matern_nuggets_stationary";
        num_params = 4;
    } else if (*kernel == 2) {
        data.kernel_fun = "bivariate_matern_flexible";
        num_params = 13;
    } else if (*kernel == 3) {
        data.kernel_fun = "bivariate_matern_parsimonious";
        num_params = 6;
    } else if (*kernel == 4) {
        data.kernel_fun = "trivariate_matern_parsimonious";
        num_params = 10;
    } else if (*kernel == 5) {
        data.kernel_fun = "univariate_spacetime_matern_stationary";
        num_params = 7;
    } else if (*kernel == 6) {
        data.kernel_fun = "bivariate_spacetime_matern_stationary";
        num_params = 10;
    } else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    //Memory allocation
    starting_theta = (double* ) malloc(num_params * sizeof(double));

    //initialize globalthetaout.
    for (j = 0; j < num_params; j++)
        *(globalthetaout + j) = -1;

    //Set data struct values based on inputs
    if (comp_mode == 0) {
        data.computation = "exact";
        data.diag_thick = 0;
    } else if (comp_mode == 1) {
#if defined( EXAGEOSTAT_USE_HICMA )
        data.computation = "lr_approx";
        data.hicma_acc = *lr_acc;
        data.hicma_maxrank = *lr_maxrank;
#endif
    } else if (comp_mode == 2) {
        data.computation = "diag_approx";
        data.diag_thick = *diag_thick;
    } else if (comp_mode == 3) {
        data.computation = "exact";
        data.diag_thick = *diag_thick;
    }

    data.dm = metric_mode == 0 ? "ed" : "gcd";


#if defined( EXAGEOSTAT_USE_HICMA )
    data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT;
#endif

    data.verbose = verbose;
    data.async = async;
    data.l1.x = x;
    data.l1.y = y;
    //data.l2            = data.l1;
    data.iter_count = 0;
    data.log = log;
    data.precision = 0;   //should be modified to be an input.

    //copy clb array to start estimation with
    for (i = 0; i < num_params; i++) {
        starting_theta[i] = clb[i];
    }

    //Optimizer initialization
    opt = nlopt_create(NLOPT_LN_BOBYQA, num_params);  //NLOPT_LN_BOBYQA  - NLOPT_GN_ORIG_DIRECT
    init_optimizer(&opt, clb, cub, *opt_tol);
    nlopt_set_maxeval(opt, *opt_max_iters);

    //Create descriptors
    if (strcmp(data.computation, "exact") == 0)
        if (data.diag_thick == 0)
            EXAGEOSTAT_dmle_Call(&data, *ncores, *gpus, *lts, *p_grid, *q_grid, *n, 0, 0);
        else {
            EXAGEOSTAT_sdmle_Call(&data, *ncores, *gpus, *lts, *p_grid, *q_grid, *n, 0, 0);
            EXAGEOSTAT_MLE_sdregister_Tile(&data);
        }
    else if (strcmp(data.computation, "diag_approx") == 0) {
        EXAGEOSTAT_dmle_diag_Call(&data, *ncores, *gpus, *lts, *p_grid, *q_grid, *n, 0, 0);
    }
#if defined( EXAGEOSTAT_USE_HICMA )
    else if (strcmp(data.computation, "lr_approx") == 0) {
        EXAGEOSTAT_TLR_dmle_Call(&data, *ncores, *gpus, *lts, *p_grid, *q_grid, *n, 0, 0);
    }
#endif

    printf("%s- %s\n", data.computation, __func__);

    double* Nrand = (double* ) malloc((*n) * 1 * sizeof(double));
    data.iter_count = 0;
    EXAGEOSTAT_TLR_MLE_dzvg_Tile(&data, Nrand, starting_theta, *n, *dts, log, *p_grid, *q_grid);

    //main algorithm call
    START_TIMING(data.total_exec_time);
    nlopt_set_max_objective(opt, MLE_alg, (void *) &data);
    nlopt_optimize(opt, starting_theta, &opt_f);
    STOP_TIMING(data.total_exec_time);

    //Destory descriptors & free memory
    nlopt_destroy(opt);

#if defined( EXAGEOSTAT_USE_HICMA )
    HICMA_Desc_Destroy((HICMA_desc_t **) &data.hicma_descC);
    HICMA_Desc_Destroy((HICMA_desc_t **) &data.hicma_descCD);
    HICMA_Desc_Destroy((HICMA_desc_t **) &data.hicma_descCUV);
    HICMA_Desc_Destroy((HICMA_desc_t **) &data.hicma_descCrk);
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descZ);
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descproduct);
    HICMA_Desc_Destroy((HICMA_desc_t **) &data.hicma_descdet);
#endif
    //copy local vector to global vector in R memory space
    for (j = 0; j < num_params; j++)
        *(globalthetaout + j) = *(starting_theta + j);

    *(globalthetaout + num_params) = data.total_exec_time / (double) data.iter_count;
    *(globalthetaout + num_params + 1) = data.total_exec_time;
    *(globalthetaout + num_params + 2) = (double) data.iter_count;
}

static void mle_general(int *kernel, double* x, int *xlen, double* y,
                        int *ylen, double* z, int *zlen,
                        double* clb, int *clblen, double* cub,
                        int *cublen, int *computation, int *diag_thick,
                        int *lr_acc, int *lr_maxrank, int *dmetric,
                        int *n, double* opt_tol, int *opt_max_iters,
                        int *ncores, int *gpus, int *ts,
                        int *p_grid, int *q_grid,
                        double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] kernel_fun:       Pointer to the stationary/nono-stationary kernel.
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector).
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] diag_thick:       Pointer to the used diagonal thick in the case of DST computation..
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] opt_tol:        Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
 * */
{
    //initialization
    int i = 0, j = 0;
    //double time_opt = 0.0, max_loglik = 0.0;
    //Log is 0 in R
    int log = 0;
    //Verbose is 0 in R
    int verbose = 0;
    //Async is 0 in R
    int async = 0;
    int comp_mode = *computation;
    int metric_mode = *dmetric;
    double opt_f;
    nlopt_opt opt;
    MLE_data data;
    double* starting_theta;

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    //Initialize the data struct
    init_data_values(&data);
    int num_params = 3;
    if (*kernel == 0) {
        data.kernel_fun = "univariate_matern_stationary";
        num_params = 3;
    } else if (*kernel == 1) {
        data.kernel_fun = "univariate_matern_nuggets_stationary";
        num_params = 4;
    } else if (*kernel == 2) {
        data.kernel_fun = "bivariate_matern_flexible";
        num_params = 13;
    } else if (*kernel == 3) {
        data.kernel_fun = "bivariate_matern_parsimonious";
        num_params = 6;
    } else if (*kernel == 4) {
        data.kernel_fun = "trivariate_matern_parsimonious";
        num_params = 10;
    } else if (*kernel == 5) {
        data.kernel_fun = "univariate_spacetime_matern_stationary";
        num_params = 7;
    } else if (*kernel == 6) {
        data.kernel_fun = "bivariate_spacetime_matern_stationary";
        num_params = 10;
    } else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    //Memory allocation
    starting_theta = (double* ) malloc(num_params * sizeof(double));

    //initialize globalthetaout.
    for (j = 0; j < num_params; j++)
        *(globalthetaout + j) = -1;

    //Set data struct values based on inputs
    if (comp_mode == 0) {
        data.computation = "exact";
        data.diag_thick = 0;
    } else if (comp_mode == 1) {
#if defined( EXAGEOSTAT_USE_HICMA )
        data.computation = "lr_approx";
        data.hicma_acc = *lr_acc;
        data.hicma_maxrank = *lr_maxrank;
#endif
    } else if (comp_mode == 2) {
        data.computation = "diag_approx";
        data.diag_thick = *diag_thick;
    } else if (comp_mode == 3) {
        data.computation = "exact";
        data.diag_thick = *diag_thick;
    }


    data.dm = metric_mode == 0 ? "ed" : "gcd";

#if defined( EXAGEOSTAT_USE_HICMA )
    data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
#endif

    data.verbose = verbose;
    data.async = async;
    data.l1.x = x;
    data.l1.y = y;
    //data.l2            = data.l1;
    data.iter_count = 0;
    data.log = log;
    data.precision = 0;   //should be modified to be an input.

    //copy clb array to start estimation with
    for (i = 0; i < num_params; i++) {
        starting_theta[i] = clb[i];
    }

    //Optimizer initialization
    opt = nlopt_create(NLOPT_LN_BOBYQA, num_params);  //NLOPT_LN_BOBYQA  - NLOPT_GN_ORIG_DIRECT
    init_optimizer(&opt, clb, cub, *opt_tol);
    nlopt_set_maxeval(opt, *opt_max_iters);

    //Create descriptors
    if (strcmp(data.computation, "exact") == 0)
        if (data.diag_thick == 0)
            EXAGEOSTAT_dmle_Call(&data, *ncores, *gpus, *ts, *p_grid, *q_grid, *n, 0, 0);
        else {
            EXAGEOSTAT_sdmle_Call(&data, *ncores, *gpus, *ts, *p_grid, *q_grid, *n, 0, 0);
            EXAGEOSTAT_MLE_sdregister_Tile(&data);
        }
    else if (strcmp(data.computation, "diag_approx") == 0) {
        EXAGEOSTAT_dmle_diag_Call(&data, *ncores, *gpus, *ts, *p_grid, *q_grid, *n, 0, 0);
    }
#if defined( EXAGEOSTAT_USE_HICMA )
    else if (strcmp(data.computation, "lr_approx") == 0) {
        EXAGEOSTAT_TLR_dmle_Call(&data, *ncores, *gpus, *ts, *p_grid, *q_grid, *n, 0, 0);
        CHAM_desc_t *CHAM_descZ = NULL;
        CHAMELEON_Desc_Create(&CHAM_descZ, NULL, ChamRealDouble, *ts, *ts, *ts * *ts, *n, 1, 0, 0, *n, 1, *p_grid,
                              *q_grid);
        data.descZ = CHAM_descZ;
    }
#endif

    printf("%s- %s\n", data.computation, __func__);
    //Copy z to descriptor
    CHAMELEON_Lapack_to_Tile(z, *n, data.descZ);
    //main algorithm call
    START_TIMING(data.total_exec_time);
    nlopt_set_max_objective(opt, MLE_alg, (void *) &data);
    nlopt_optimize(opt, starting_theta, &opt_f);
    STOP_TIMING(data.total_exec_time);

    //Destory descriptors & free memory
    nlopt_destroy(opt);
    if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "diag_approx") == 0) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descC);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descZ);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descZcpy);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descproduct);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.descdet);
    }
#if defined( EXAGEOSTAT_USE_HICMA )
    else if (strcmp(data.computation, "tlr_approx") == 0) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descCD);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descCUV);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descCrk);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descZ);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descZcpy);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descproduct);
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &data.hicma_descdet);

    }
#endif
    //copy local vector to global vector in R memory space
    for (j = 0; j < num_params; j++)
        *(globalthetaout + j) = *(starting_theta + j);

    *(globalthetaout + num_params) = data.total_exec_time / (double) data.iter_count;
    *(globalthetaout + num_params + 1) = data.total_exec_time;
    *(globalthetaout + num_params + 2) = (double) data.iter_count;
}


void mle_exact(double* x, int *xlen, double* y,
               int *ylen, double* z, int *zlen,
               double* clb, int *clblen, double* cub,
               int *cublen, int *kernel, int *dmetric,
               int *n, double* opt_tol, int *opt_max_iters,
               int *ncores, int *gpus, int *ts,
               int *p_grid, int *q_grid, double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector)
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] kernel:           Pointer to the computation kernel.
 * @param[in] kernellen:        Pointer to the length of computation kernel.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3)
 * */
{
    int comp_mode = 0;  //exact

    mle_general(kernel, x, xlen, y, ylen,
                z, zlen, clb, clblen,
                cub, cublen, &comp_mode,
                0, 0, 0,
                dmetric, n, opt_tol, opt_max_iters,
                ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}


void mle_exact_non_stat(double* x, int *xlen, double* y,
                        int *ylen, double* z, int *zlen,
                        int *kernel, double* clb, int *clblen,
                        double* cub, int *cublen, int *dmetric,
                        int *n, double* opt_tol, int *opt_max_iters,
                        int *ncores, int *gpus, int *ts,
                        int *p_grid, int *q_grid, double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] kernel:           Pointer to the stationary/nono-stationary kernel.
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector)
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] kernel:           Pointer to the computation kernel.
 * @param[in] kernellen:        Pointer to the length of computation kernel.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3)
 * */
{
    int comp_mode = 0;
    mle_general(kernel, x, xlen, y, ylen,
                z, zlen, clb, clblen,
                cub, cublen, &comp_mode,
                0, 0, 0,
                dmetric, n, opt_tol, opt_max_iters,
                ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}

void mle_tlr(double* x, int *xlen, double* y,
             int *ylen, double* z, int *zlen,
             int *kernel, double* clb, int *clblen,
             double* cub, int *cublen, int *tlr_acc,
             int *tlr_maxrank, int *dmetric, int *n,
             double* opt_tol, int *opt_max_iters, int *ncores,
             int *gpus,  int *dts, int *lts, int *p_grid,
             int *q_grid, double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector)
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector).
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] kernel:           Pointer to the computation kernel.
 * @param[in] kernellen:        Pointer to the length of computation kernel.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
 * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
 * */
{
    int comp_mode = 1;

    mle_general_tlr(kernel, x, xlen, y, ylen, z, zlen,
                    clb, clblen, cub, cublen,
                    &comp_mode, 0, tlr_acc, tlr_maxrank,
                    dmetric, n, opt_tol, opt_max_iters,
                    ncores, gpus, dts, lts, p_grid, q_grid, globalthetaout);
}


void mle_dst(double* x, int *xlen, double* y,
             int *ylen, double* z, int *zlen,
             int *kernel, double* clb, int *clblen,
             double* cub, int *cublen, int *dst_thick,
             int *dmetric, int *n, double* opt_tol,
             int *opt_max_iters, int *ncores, int *gpus,
             int *ts, int *p_grid, int *q_grid,
             double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector).
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] kernel:           Pointer to the computation kernel.
 * @param[in] kernellen:        Pointer to the length of computation kernel.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
 * */
{
    int comp_mode = 2;

    mle_general(kernel, x, xlen, y, ylen,
                z, zlen, clb, clblen,
                cub, cublen, &comp_mode,
                dst_thick, 0, 0,
                dmetric, n, opt_tol, opt_max_iters,
                ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}


void mle_mp(double* x, int *xlen, double* y,
            int *ylen, double* z, int *zlen,
            int *kernel, double* clb, int *clblen,
            double* cub, int *cublen, int *mp_band,
            int *dmetric, int *n, double* opt_tol,
            int *opt_max_iters, int *ncores, int *gpus,
            int *ts, int *p_grid, int *q_grid,
            double* globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] x:                Pointer to the x vector.
 * @param[in] xlen:             Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector).
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] kernel:           Pointer to the computation kernel.
 * @param[in] kernellen:        Pointer to the length of computation kernel.
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
 * */
{
    int comp_mode = 3;

    mle_general(kernel, x, xlen, y, ylen,
                z, zlen, clb, clblen,
                cub, cublen, &comp_mode,
                mp_band, 0, 0,
                dmetric, n, opt_tol, opt_max_iters,
                ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}


void exact_predict(double* xobs, int *xobs_len, double* yobs,
                   int *yobs_len, double* zobs, int *zobs_len,
                   double* xmiss, int *xmiss_len, double* ymiss,
                   int *ymiss_len, int *nZobs, int *nZmiss,
                   int *kernel, double* est_theta, int *thetalen,
                   int *computation, int *dmetric, int *ncores,
                   int *gpus, int *ts, int *p_grid,
                   int *q_grid, double* globalzpredict)
//! R-wrapper to predict missing values on certain missing locations.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector
 * @param[in] xobs:             Pointer to the xobs vector.
 * @param[in] xobs_len:         Pointer to the length of xobs vector.
 * @param[in] yobs:             Pointer to the yobs vector.
 * @param[in] yobs_len:         Pointer to the length of yobs vector.
 * @param[in] z:                Pointer to the zobs vector (measurements).
 * @param[in] nZobs:         Pointer to the length of zobs vector (measurements).
 * @param[in] xmiss:            Pointer to the xmiss vector.
 * @param[in] xmiss_len:        Pointer to the length of xmiss vector.
 * @param[in] ymiss:            Pointer to the ymiss vector.
 * @param[in] ymiss_len:        Pointer to the length of ymiss vector.
 * @param[in] z:                Pointer to the zmiss vector (measurements).
 * @param[in] nZmiss:        Pointer to the length of zmiss vector (measurements).
 * @param[in] theta:            Pointer to the theta vector.
 * @param[in] thetalen:         Pointer to the length of y-dim vector.
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
 * */
{
    //initialization
    int i = 0, j = 0, log = 0, verbose = 0, async = 0, num_params = 0;
    double opt_f;
    nlopt_opt opt;
    MLE_data data;
    double* starting_theta, *theta;
    int p = 1;                        //number of variables.
    double* zmiss;
    char *kernel_fun;
    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }
    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    //Initialize the data struct
    init_data_values(&data);

    int comp_mode = *computation;
    int metric_mode = *dmetric;

    theta = (double* ) malloc(num_params * sizeof(double));

    //Assign initial_ theta vector to generate synthetic dataset
    for (j = 0; j < *thetalen; j++)
        theta[j] = est_theta[j];

    num_params = get_num_params(data.kernel_fun);

    CHAMELEON_Sequence_Create(&msequence);
    data.sequence = msequence;
    data.request = mrequest;

    //initialize globalthetaout.
    for (j = 0; j < (*nZmiss); j++)
        *(globalzpredict + j) = -1000;

    //Set data struct values based on inputs
    if (comp_mode == 0) {
        data.computation = "exact";
        data.diag_thick = 0;
    } else
        printf("In R mode: only exact prediction is allowed!\n");

    data.dm = metric_mode == 0 ? "ed" : "gcd";

#if defined( EXAGEOSTAT_USE_HICMA )
    data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
#endif

    data.verbose = verbose;
    data.async = async;
    data.lobs.x = xobs;
    data.lobs.y = yobs;
    data.lmiss.x = xmiss;
    data.lmiss.y = ymiss;
    data.log = log;

    if (*kernel == 0)
        data.kernel_fun = "univariate_matern_stationary";
    else if (*kernel == 1)
        data.kernel_fun = "univariate_matern_nuggets_stationary";
    else if (*kernel == 2)
        data.kernel_fun = "bivariate_matern_flexible";
    else if (*kernel == 3)
        data.kernel_fun = "bivariate_matern_parsimonious";
    else if (*kernel == 4)
        data.kernel_fun = "trivariate_matern_parsimonious";
    else if (*kernel == 5)
        data.kernel_fun = "univariate_spacetime_matern_stationary";
    else if (*kernel == 6)
        data.kernel_fun = "bivariate_spacetime_matern_stationary";
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    data.precision = 0;               //should be modified to be an input.
    zmiss = (double* ) calloc(p * (*nZmiss), sizeof(double));

    if (strcmp(data.computation, "exact") == 0)
        prediction_init(&data, (*nZmiss), (*nZobs), (*ts), (*p_grid), (*q_grid), 0);

    CHAM_desc_t *CHAM_descZobs = (CHAM_desc_t *) (data.descZobs);
    CHAM_desc_t *CHAM_descZmiss = (CHAM_desc_t *) (data.descZmiss);
    CHAMELEON_Lapack_to_Tile(zobs, (*nZobs), CHAM_descZobs);

    if (strcmp(data.computation, "exact") == 0)
        EXAGEOSTAT_dmle_Predict_Tile(&data, theta, p * (*nZmiss),
                                     p * (*nZobs), zobs, NULL, zmiss, (*nZobs));

    CHAMELEON_Tile_to_Lapack(CHAM_descZmiss, zmiss, (*nZmiss));

    //copy local vector to global vector in R memory space
    for (j = 0; j < (*nZmiss); j++)
        *(globalzpredict + j) = zmiss[j];
}

void exact_mloe_mmom(double* xobs, int *xobs_len, double* yobs,
                     int *yobs_len, double* zobs, int *zobs_len,
                     double* xmiss, int *xmiss_len, double* ymiss,
                     int *ymiss_len, int *nZobs, int *nZmiss,
                     int *kernel, double* est_theta, double* true_theta,
                     int *thetalen, int *dmetric,
                     int *ncores, int *gpus, int *ts,
                     int *p_grid, int *q_grid, double* global_mloe_mmom) {

    //initialization
    int i = 0, j = 0, log = 0, verbose = 0, async = 0, num_params = 0;
    double opt_f;
    nlopt_opt opt;
    MLE_data data;
    double* starting_theta;
    int p = 1;                        //number of variables.
    double* zmiss;
    char *kernel_fun;
    CHAM_context_t *CHAM;
    CHAM = chameleon_context_self();
    if (CHAM == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }
    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

    //Initialize the data struct
    init_data_values(&data);

    int comp_mode = 0;    //exact
    int metric_mode = *dmetric;
    double* theta = (double* ) malloc(num_params * sizeof(double));


    CHAMELEON_Sequence_Create(&msequence);
    data.sequence = msequence;
    data.request = mrequest;

    //initialize globalthetaout.
    for (j = 0; j < 2; j++)
        *(global_mloe_mmom + j) = -1000;

    //Set data struct values based on inputs
    if (comp_mode == 0) {
        data.computation = "exact";
        data.diag_thick = 0;
    } else
        printf("In R mode: only exact prediction is allowed!\n");

    data.dm = metric_mode == 0 ? "ed" : "gcd";

#if defined( EXAGEOSTAT_USE_HICMA )
    data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
#endif

    data.verbose = verbose;
    data.async = async;
    data.lobs.x = xobs;
    data.lobs.y = yobs;
    data.lmiss.x = xmiss;
    data.lmiss.y = ymiss;
    data.log = log;

    if (*kernel == 0)
        data.kernel_fun = "univariate_matern_stationary";
    else if (*kernel == 1)
        data.kernel_fun = "univariate_matern_nuggets_stationary";
    else if (*kernel == 2)
        data.kernel_fun = "bivariate_matern_flexible";
    else if (*kernel == 3)
        data.kernel_fun = "bivariate_matern_parsimonious";
    else if (*kernel == 4)
        data.kernel_fun = "trivariate_matern_parsimonious";
    else if (*kernel == 5)
        data.kernel_fun = "univariate_spacetime_matern_stationary";
    else if (*kernel == 6)
        data.kernel_fun = "bivariate_spacetime_matern_stationary";
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    data.precision = 0;               //should be modified to be an input.
    zmiss = (double* ) calloc(p * (*nZmiss), sizeof(double));


    mloe_mmom_init(&data, (*nZmiss), (*nZobs), (*ts), (*p_grid), (*q_grid));

    EXAGEOSTAT_dmle_mloe_mmom_Tile(&data, true_theta, est_theta,
                                   p * (*nZmiss), p * (*nZobs), (*nZobs + *nZmiss));

    //copy local vector to global vector in R memory space
    *(global_mloe_mmom) = data.mloe;
    *(global_mloe_mmom + 1) = data.mmom;
}


void rexageostat_init(int *ncores, int *gpus, int *dts, int *lts)
//! R-wrapper to initiate exageostat.
/*!  -- using dense or approximate computation
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not CHAM.
 * */
{
    exageostat_init(ncores, gpus, dts, lts);  ///this call shoudl be modified to include both dts and lts.
}

void rexageostat_finalize()
//! R-wrapper to finalize exageostat.
{
    exageostat_finalize();
}


void fisher_general(double* x, double* y, int *n,
                    double* est_theta, int *thetalen,
                    int *ts, int *dmetric, int *p_grid, int *q_grid,
                    double* globalthetaout) {

    MLE_data data;
    int j = 0;
    double* theta;
    RUNTIME_sequence_t *msequence;
    int metric_mode = *dmetric;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&msequence);

    //Initialize the data struct
    init_data_values(&data);

    data.l1.x = x;
    data.l1.y = y;
    data.sequence = msequence;
    data.request = mrequest;
    data.dm = metric_mode == 0 ? "ed" : "gcd";

    theta = (double* ) malloc(*thetalen * sizeof(double));

    for (j = 0; j < *thetalen; j++)
        theta[j] = est_theta[j];

    data.ooc = 0;

}