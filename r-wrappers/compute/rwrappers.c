/**
 *
 * Copyright (c) 2017-2019  King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2019-07-30
 *
 **/
#include "../include/rwrappers.h"
void  gen_z_givenlocs_exact(double *x, int *xlen, double *y,
        int *ylen,  double *theta1,
        double *theta2, double *theta3,
        int *dmetric, int *n,
        int *ncores, int *gpus, int *ts,
        int *p_grid, int *q_grid, int *veclen,
        double *globalvec)
    //! direct function to generate synthetics datasets (X, Y) 2D locations Z measurement vector.
    /*!  -- using dense or approximate computation
     * Returns Z observation vector
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] x:                Pointer to the x vector.
     * @param[in] xlen:             Pointer to the length of x vector.
     * @param[in] y:                Pointer to the y vector.
     * @param[in] ylen:             Pointer to the length of y vector.
     * @param[in] theta1:           Pointer to the variance value (theta1).
     * @param[in] theta2:           Pointer to the range value (theta2).
     * @param[in] theta3:           Pointer to the smoothness value (theta3).
     * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx).
     * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
     * @param[in] n:                Pointer to the problem size (number spatial locations).
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] veclen:           Pointer to the length of global vec (R memory space).
     * @param[in] globalvec:        Pointer to R memory space (x:y:z).
     * */
{
    //Initialization
    //int i            = 0;
    int j               = 0;
    //Log is 0 in R
    int log                = 0;
    //Verbose is 0 in R
    int verbose         = 0;
    //Async is 0 in R
    int async           = 0;
    //exact -> 0 (Only exact is needed)
    int comp_mode       = 0;
    int metric_mode     = 0;
    double  *initial_theta, *localvec;
    MLE_data data;
    location *locations;
    int iseed[4]={0, 0, 0, 1};
    MORSE_sequence_t *msequence;
    MORSE_request_t mrequest[2] = { MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER };

    MORSE_context_t *morse;
    morse = morse_context_self();
    if (morse == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    //initialize globalvec
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = 0;    
    MORSE_desc_t *MORSE_descC       = NULL;
    MORSE_desc_t *MORSE_descZ       = NULL;

    //Memory allocation
    locations         = (location *) malloc( sizeof(location*));
    initial_theta           = (double *) malloc(3 * sizeof(double));
    localvec                = (double *) malloc( *veclen * sizeof(double));

    //Assign x, y vectors to local space (C memeory space)
    //data.l1.x             = localvec;
    //data.l1.y             = &localvec[*n];

    //set XY locs
    locations->x        = x;
    locations->y        = y;
    data.l1                 = *locations;

    //Assign x, y vectors to local space (C memeory space)
    /*for(i = 0; i < *n; i++)
      localvec[i] = data.l1.x[i];
      for(i = *n; i < 2*(*n); i++)
      localvec[i] = data.l1.y[i-(*n)];
      */

    //Set data struct values based on inputs
    comp_mode               = 0;  //Usually exact
    data.computation    = comp_mode == 0 ? "exact" : "appro";
    data.dm                 = metric_mode == 0 ? "ed" : "gcd";
    data.verbose            = verbose;
    //data.l2                       = data.l1;
    data.async              = async;
    //data.log                = log;
    data.obsFPath           = "";
    data.log                = log;
    gsl_set_error_handler_off () ;

    //Assign initial_ theta vector to generate synthetic dataset
    initial_theta[0]        = *theta1;
    initial_theta[1]        = *theta2;
    initial_theta[2]        = *theta3;
    //Nomral random generation of e -- ei~N(0, 1) to generate Z
    LAPACKE_dlarnv(3, iseed, *n, &localvec[0]);


    //Create Descriptors
    MORSE_Sequence_Create(&msequence);
    MORSE_Desc_Create(&MORSE_descC, NULL , MorseRealDouble, *ts, *ts, *ts * *ts, *n, *n, 0, 0, *n, *n, *p_grid, *q_grid);
    MORSE_Desc_Create(&MORSE_descZ, NULL, MorseRealDouble, *ts, *ts, *ts * *ts, *n, 1,  0, 0, *n , 1, *p_grid, *q_grid);
    data.sequence          = msequence;
    data.request           = mrequest;

    data.descC       = MORSE_descC;
    data.descZ       = MORSE_descZ;

    //Main algorithm call
    MLE_zvg(&data, localvec, initial_theta, *n, *ts, log, *p_grid, *q_grid) ;

    MORSE_Tile_to_Lapack(data.descZ, localvec, *n);
    //copy local vector to global vector in R memory space
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = *(localvec + j);

    //Destory descriptors
    MORSE_Desc_Destroy((MORSE_desc_t **) &data.descC );
    MORSE_Desc_Destroy((MORSE_desc_t **) &data.descZ);

    //free memory
    free(initial_theta);
    free(localvec);
}
void  gen_z_exact( double *theta1, double *theta2,
        double *theta3, int *dmetric, int *n,
        int *seed,
        int *ncores,  int *gpus,
        int *ts, int *p_grid, int *q_grid,
        int *veclen,  double *globalvec)
    //! direct function to generate synthetics datasets (X, Y) 2D locations Z measurement vector.
    /*!  -- using dense or approximate computation
     * Returns Z observation vector
     * @param[in] n:                Pointer to the problem size (number spatial locations).
     * @param[in] theta1:           Pointer to the variance value (theta1).
     * @param[in] theta2:           Pointer to the range value (theta2).
     * @param[in] theta3:           Pointer to the smoothness value (theta3).
     * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
     * @param[in] seed:            Pointer to the random seed.    
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] veclen:           Pointer to the length of global vec (R memory space).
     * @param[in] globalvec:        Pointer to R memory space (x:y:z).
     * */
{
    //Initialization
    int i            = 0;
    int j               = 0;
    //Log is 0 in R
    int log                = 0;
    //Verbose is 0 in R
    int verbose         = 0;
    //Async is 0 in R
    int async           = 0;
    //exact -> 0 (Only exact is needed)
    int comp_mode       = 0;
    int metric_mode     = 0;
    double  *initial_theta, *localvec;
    MLE_data data;
    location *locations;
    int iseed[4]={*seed, *seed, *seed, 1};
    MORSE_sequence_t *msequence;
    MORSE_request_t mrequest[2] = { MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER };

    MORSE_context_t *morse;
    morse = morse_context_self();
    if (morse == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    //initialize globalvec
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = 0;

    MORSE_desc_t *MORSE_descC       = NULL;
    MORSE_desc_t *MORSE_descZ       = NULL;

    //Memory allocation
    initial_theta           = (double *) malloc(3 * sizeof(double));
    localvec                = (double *) malloc( *veclen * sizeof(double));

    //Assign x, y vectors to local space (C memeory space)
    //data.l1.x             = localvec;
    //data.l1.y             = &localvec[*n];

    //Generate XY locs
    locations               = GenerateXYLoc(*n, *seed);
    //printf("####seed..........: %d\n", *seed);
    data.l1                 = *locations;

    //Assign x, y vectors to local space (C memeory space)
    for(i = 0; i < *n; i++)
        localvec[i] = data.l1.x[i];
    for(i = *n; i < 2*(*n); i++)
        localvec[i] = data.l1.y[i-(*n)];

    //Set data struct values based on inputs
    comp_mode               = 0;  //Usually exact
    data.computation    = comp_mode == 0 ? "exact" : "appro";
    data.dm                 = metric_mode == 0 ? "ed" : "gcd";
    data.verbose            = verbose;
    //data.l2                       = data.l1;
    data.async              = async;
    //data.log                = log;
    data.obsFPath           = "";
    data.log                = log;
    data.precision           = 0;   
    gsl_set_error_handler_off () ;

    //Assign initial_ theta vector to generate synthetic dataset
    initial_theta[0]        = *theta1;
    initial_theta[1]        = *theta2;
    initial_theta[2]        = *theta3;
    //Nomral random generation of e -- ei~N(0, 1) to generate Z
    LAPACKE_dlarnv(3, iseed, *n, &localvec[2**n]);

    //Create Descriptors
    MORSE_Sequence_Create(&msequence);
    MORSE_Desc_Create(&MORSE_descC, NULL , MorseRealDouble, *ts, *ts, *ts * *ts, *n, *n, 0, 0, *n, *n, *p_grid, *q_grid);
    MORSE_Desc_Create(&MORSE_descZ, NULL, MorseRealDouble, *ts, *ts, *ts * *ts, *n, 1,  0, 0, *n , 1, *p_grid, *q_grid);
    data.sequence          = msequence;
    data.request           = mrequest;


    data.descC    = MORSE_descC;
    data.descZ    = MORSE_descZ;

    //Main algorithm call
    MLE_zvg(&data, &localvec[2**n], initial_theta, *n, *ts, log, *p_grid, *q_grid) ;

    MORSE_Tile_to_Lapack(data.descZ, &localvec[2**n], *n);
    //copy local vector to global vector in R memory space
    for (j = 0; j < *veclen; j++)
        *(globalvec + j) = *(localvec + j);

    //Destory descriptors
    MORSE_Desc_Destroy((MORSE_desc_t **) &data.descC );
    MORSE_Desc_Destroy((MORSE_desc_t **) &data.descZ);

    //free memory
    free(initial_theta);
    free(localvec);
}

static void  mle_general(char *kernel_fun, double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen,  int *computation,  int *diag_thick,
        int *lr_acc, int *lr_maxrank,  int *dmetric,
        int *n,  double *opt_tol, int *opt_max_iters,
        int *ncores, int *gpus,    int *ts,
        int *p_grid, int *q_grid,
        double *globalthetaout)
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
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
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
    int log                = 0;
    //Verbose is 0 in R
    int verbose         = 0;        
    //Async is 0 in R
    int async           = 0;
    int comp_mode       = *computation;
    int metric_mode     = *dmetric;
    double opt_f;
    nlopt_opt opt;
    MLE_data data;
    double *starting_theta;

    MORSE_context_t *morse;
    morse = morse_context_self();
    if (morse == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
        return;
    }

    init_data_values(&data);

    int  num_params;
    if(strcmp(data.kernel_fun, "univariate_matern_stationary")   == 0)
        num_params = 3;
    else if(strcmp(data.kernel_fun, "univariate_matern_non_stationary")   == 0)
        num_params = 9;
    else if(strcmp(data.kernel_fun, "bivariate_matern_flexible")   == 0)
        num_params = 13;
    else if(strcmp(data.kernel_fun, "bivariate_matern_parsimonious")   == 0)
        num_params = 6;
    else if(strcmp(data.kernel_fun, "bivariate_matern_parsimonious2")   == 0)
        num_params = 6;
    else if(strcmp(data.kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
        num_params = 6;
    else if(strcmp(data.kernel_fun, "univariate_spacetime_matern_stationary")   == 0)
        num_params = 7;
    else
    {
        fprintf(stderr,"Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n",__func__);
        exit(0);
    }
    //num_params = strcmp(data->kernel_fun, "stationary_kernel")   == 0? 3 : 9;

    //Memory allocation
    starting_theta         = (double *) malloc(num_params * sizeof(double));

    //initialize globalthetaout.
    for (j = 0; j < num_params; j++)
        *(globalthetaout + j) = -1;

    //Set data struct values based on inputs
    if(comp_mode == 0)
        data.computation        = "exact";
    else if(comp_mode == 1)
    {
#if defined( EXAGEOSTAT_USE_HICMA )
        data.computation        = "lr_approx";
        data.hicma_acc          = *lr_acc;
        data.hicma_maxrank      = *lr_maxrank;
#endif
    }
    else if(comp_mode == 2)
    {
        data.computation        = "diag_approx";
        data.diag_thick         = *diag_thick;
    }

    data.dm         = metric_mode == 0? "ed" : "gcd";
    //if (strcmp (data.dm, "ed") == 0)
    //        data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT;
    //else

#if defined( EXAGEOSTAT_USE_HICMA )
    data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;    
#endif

    data.verbose        = verbose;
    data.async          = async;
    data.l1.x           = x;
    data.l1.y           = y;
    //data.l2            = data.l1;    
    data.iter_count     = 0;
    data.log            = log;
    data.kernel_fun     = kernel_fun;
    data.precision      = 0;   //should be modified to be an input.

    //copy clb array to start estimation with
    for(i=0;i<num_params;i++)
        starting_theta[i] = clb[i];

    //Optimizer initialization
    opt=nlopt_create( NLOPT_LN_BOBYQA, num_params);  //NLOPT_LN_BOBYQA  - NLOPT_GN_ORIG_DIRECT
    init_optimizer(&opt, clb, cub, *opt_tol);
    nlopt_set_maxeval(opt, *opt_max_iters);


    //Create descriptors
    if(strcmp (data.computation, "exact") == 0)        
        MORSE_dmle_Call(&data, *ncores,*gpus, *ts, *p_grid, *q_grid, *n,  0, 0);
    else if (strcmp (data.computation, "diag_approx") == 0)
        MORSE_dmle_diag_Call(&data, *ncores, *gpus, *ts, *p_grid, *q_grid, *n, 0, 0);
#if defined( EXAGEOSTAT_USE_HICMA )
    else if (strcmp (data.computation, "lr_approx") == 0)
    {
#if defined( EXAGEOSTAT_USE_HICMA )
        HICMA_dmle_Call(&data, *ncores, *gpus, *ts, *p_grid, *q_grid, *n,  0, 0);
        MORSE_desc_t *MORSE_descZ       = NULL;
        MORSE_Desc_Create(&MORSE_descZ, NULL, MorseRealDouble, *ts, *ts, *ts * *ts, *n, 1, 0, 0, *n, 1, *p_grid, *q_grid);
        data.descZ = MORSE_descZ;
#endif
    }
#endif

    printf("%s- %s\n", data.computation, __func__);

    //Copy z to descriptor
    MORSE_Lapack_to_Tile( z, *n, data.descZ);

    //print summary message
    //print_summary(1, *n, *ncores, *gpus, *ts, data.computation, 1, *p_grid, *q_grid);

    //main algorithm call
    START_TIMING(data.total_exec_time);
    nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
    nlopt_optimize(opt, starting_theta, &opt_f);
    STOP_TIMING(data.total_exec_time);

    //Print results to log files
    //print_result(&data, starting_theta, *n, 1, *ncores, *ts, 1, NULL, computation, 1, 1, data.final_loglik);

    //Destory descriptors & free memory
    nlopt_destroy(opt);



    if(strcmp (data.computation, "exact") == 0 || strcmp (data.computation, "diag_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.descC );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.descZ );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.descZcpy );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.descproduct );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.descdet );
    }
#if defined( EXAGEOSTAT_USE_HICMA )
    else if(strcmp (data.computation, "tlr_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descCD );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descCUV );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descCrk );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descZ );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descZcpy );        
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descproduct );
        MORSE_Desc_Destroy( (MORSE_desc_t**) &data.hicma_descdet );        

    }
#endif
    //copy local vector to global vector in R memory space
    for (j = 0; j < num_params; j++)
        *(globalthetaout + j) = *(starting_theta + j);

    *(globalthetaout + num_params) = data.total_exec_time/(double)data.iter_count;
    *(globalthetaout + num_params+1) = data.total_exec_time; 
    *(globalthetaout + num_params+2) = (double)data.iter_count;
}


void  mle_exact(double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *dmetric, int *n, double *opt_tol,
        int *opt_max_iters, int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout)
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
     * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
     * @param[in] n:                Pointer to the problem size (number spatial locations).
     * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
     * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3)
     * */
{
    int comp_mode = 0;

    mle_general( "univariate_matern_stationary", x, xlen, y, ylen,
            z, zlen, clb, clblen, 
            cub, cublen, &comp_mode,
            0, 0, 0,
            dmetric, n, opt_tol, opt_max_iters,
            ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}


void  mle_exact_non_stat(double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *dmetric, int *n, double *opt_tol,
        int *opt_max_iters, int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout)
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
     * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
     * @param[in] n:                Pointer to the problem size (number spatial locations).
     * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
     * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3)
     * */
{
    int comp_mode = 0;

    mle_general( "univariate_matern_non_stationary", x, xlen, y, ylen,
            z, zlen, clb, clblen,
            cub, cublen, &comp_mode,
            0, 0, 0,
            dmetric, n, opt_tol, opt_max_iters,
            ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}

void  mle_tlr(  double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *tlr_acc, int *tlr_maxrank,
        int *dmetric, int *n,  double *opt_tol,
        int *opt_max_iters, int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout)
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
     * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
     * @param[in] n:                Pointer to the problem size (number spatial locations).
     * @param[in] opt_tol:          Pointer to the tol parameter ( tolerance that is used for the purpose of stopping criteria only).
     * @param[in] opt_max_iters:    Pointer to the maximum number of mle iterations.
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
     * */
{
    int comp_mode = 1;
    mle_general( "univariate_matern_stationary", x, xlen, y, ylen, z, zlen,
            clb, clblen, cub, cublen,
            &comp_mode,  0, tlr_acc, tlr_maxrank,
            dmetric, n, opt_tol, opt_max_iters,
            ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}


void  mle_dst(  double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *dst_thick, int *dmetric,
        int *n, double *opt_tol, int *opt_max_iters,
        int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout)
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
     * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd).
     * @param[in] n:                Pointer to the problem size (number spatial locations).
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
     * @param[in] globalthetaout:   Pointer to R memory space (theta1:theta2:theta3).
     * */
{
    int comp_mode = 2;

    mle_general( "univariate_matern_stationary", x, xlen, y, ylen,
            z, zlen, clb, clblen,
            cub, cublen, &comp_mode,
            dst_thick, 0, 0,
            dmetric, n, opt_tol, opt_max_iters,
            ncores, gpus, ts, p_grid, q_grid, globalthetaout);
}





void rexageostat_init(int *ncores, int *gpus, int *ts)
    //! R-wrapper to initiate exageostat.
    /*!  -- using dense or approximate computation
     * @param[in] ncores:           Pointer to the number of CPUs.
     * @param[in] gpus:             Pointer to the number of GPUs.
     * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
     * */
{
    exageostat_init(ncores, gpus, ts, 0);  ///this call shoudl be modified to include both dts and lts.
}
void rexageostat_finalize()
    //! R-wrapper to finalize exageostat.
{
    exageostat_finalize();
}
