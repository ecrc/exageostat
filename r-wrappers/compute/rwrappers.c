/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-14
 *
 **/
#include "../include/rwrappers.h"

void  rexageostat_gen_z(int *n, int *ncores,  int *gpus,  int *ts,  int *p_grid, int *q_grid,  double *theta1, double *theta2, double *theta3,  int *computation, int *dmetric, int *veclen,  double *globalvec)
//! R-wrapper to generate synthetics datasets (X, Y) 2D locations Z measurement vector.
/*!  -- using dense or approximate computation
 * Returns Z observation vector
 * @param[in] n:		Pointer to the problem size (number spatial locations).
 * @param[in] ncores:		Pointer to the number of CPUs.
 * @param[in] gpus:		Pointer to the number of GPUs.
 * @param[in] ts:		Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE. 
 * @param[in] p_grid:		Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:		Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] theta1:		Pointer to the variance value (theta1).
 * @param[in] theta2:		Pointer to the range value (theta2).
 * @param[in] theta3:		Pointer to the smoothness value (theta3).
 * @param[in] computation:	Pointer to the computation mode (0--> exact, 1-->approx).
 * @param[in] dmetric:		Pointer to the used metric (0-->ed, 1-->gcd).
 * @param[in] veclen:		Pointer to the length of global vec (R memory space).
 * @param[in] globalvec:	Pointer to R memory space (x:y:z).
 * */
{
	//Initialization
        int i		= 0;
	int j		= 0;
	//Log is 0 in R
	int log		= 0;
	//Verbose is 0 in R
	int verbose	= 0;
	//Async is 0 in R
	int async	= 0;
	int comp_mode   = *computation;
	int metric_mode = *dmetric;
	double  *initial_theta, *localvec;
        MLE_data data;
        int iseed[4]={0, 0, 0, 1};
        MORSE_sequence_t *msequence;
        MORSE_request_t mrequest[2] = { MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER };

	//Memory allocation
	initial_theta		= (double *) malloc(3 * sizeof(double));
	localvec		= (double *) malloc( *veclen * sizeof(double));

	//Assign x, y vectors to local space (C memeory space)
	data.l1.x		= localvec;
	data.l1.y		= &localvec[*n];

	//Generate XY locs
	GenerateXYLoc(*n, "", &data.l1);
	
        //Set data struct values based on inputs
	data.computation	= comp_mode == 0 ? "exact" : "appro";  
        data.dm 		= comp_mode == 0 ? "ed" : "gcd";
        data.verbose		= verbose;
        data.l2			= data.l1;
	data.async		= async;
        data.log		= log;
	data.obsFPath		= "";
	data.log		= log;
        gsl_set_error_handler_off () ;
	
  	//Assign initial_ theta vector to generate synthetic dataset
	initial_theta[0]	= *theta1;
	initial_theta[1]	= *theta2;
	initial_theta[2]	= *theta3;      
	//Nomral random generation of e -- ei~N(0, 1) to generate Z
	LAPACKE_dlarnv(3, iseed, *n, &localvec[2**n]);

	//Create Descriptors
	MORSE_Sequence_Create(&msequence);
        MORSE_Desc_Create(&data.descC, NULL , MorseRealDouble, *ts, *ts, *ts * *ts, *n, *n, 0, 0, *n, *n, *p_grid, *q_grid);
        MORSE_Desc_Create(&data.descZ, NULL, MorseRealDouble, *ts, *ts, *ts * *ts, *n, 1,  0, 0, *n , 1, *p_grid, *q_grid);
	data.sequence          = msequence;
        data.request           = mrequest;

        //Main algorithm call
        MLE_zvg(&data, &localvec[2**n], initial_theta, *n, ts, 1, log) ;
	
	MORSE_Tile_to_Lapack(data.descZ, &localvec[2**n], *n);    
	//copy local vector to global vector in R memory space
	for (j = 0; j < *veclen; j++)
		*(globalvec + j) = *(localvec + j);	

	//Destory descriptors
	MORSE_Desc_Destroy(&data.descC );
        MORSE_Desc_Destroy(&data.descZ);
	
	//free memory
	free(initial_theta);
	free(localvec);
}

void  rexageostat_likelihood(int *n,  int *ncores, int *gpus, int *ts, int *p_grid, int *q_grid,  double *x, int *xlen, double *y, int *ylen, double *z, int *zlen, double *clb, int *clblen, double *cub, int *cublen,  int *computation, int *dmetric, double *globalthetaout)
//! R-wrapper to estimate the makimum likelihood function.
/*!  -- using dense or approximate computation
 * Returns the optimized theta vector 
 * @param[in] n:                Pointer to the problem size (number spatial locations).
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
 * @param[in] p_grid:           Pointer to the p_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] q_grid:           Pointer to the q_grid value ( >1 in the case of distributed systems, 1--> others).
 * @param[in] x:           	Pointer to the x vector.
 * @param[in] xlen:           	Pointer to the length of x vector.
 * @param[in] y:                Pointer to the y vector.
 * @param[in] ylen:             Pointer to the length of y vector.
 * @param[in] z:                Pointer to the z vector (measurements).
 * @param[in] zlen:             Pointer to the length of z vector (measurements).
 * @param[in] clb:              Pointer to the clb vector (lower bound optimization vector).
 * @param[in] clblen:           Pointer to the length of z vector (lower bound optimization vector).
 * @param[in] cub:              Pointer to the z vector (upper bound optimization vector) 
 * @param[in] cublen:           Pointer to the length of z vector (upper bound optimization vector).
 * @param[in] computation:      Pointer to the computation mode (0--> exact, 1-->approx)
 * @param[in] dmetric:          Pointer to the used metric (0-->ed, 1-->gcd)
 * @param[in] globalthetaout:        Pointer to R memory space (theta1:theta2:theta3)
 * */
{
	//initialization
        int i = 0, j = 0;
	double time_opt = 0.0, max_loglik = 0.0;
        //Log is 0 in R
        int log         = 0;
        //Verbose is 0 in R
        int verbose     = 0;        
        //Async is 0 in R
        int async       = 0;
	int comp_mode   = *computation;
        int metric_mode = *dmetric;
	double * opt_f;
        nlopt_opt opt;
        MLE_data data;
	double *starting_theta;


	//Memory allocation
	starting_theta		 = (double *) malloc(3 * sizeof(double));

        //Set data struct values based on inputs
        data.computation	= comp_mode == 0? "exact" : "appro";
        data.dm 		= comp_mode == 0? "ed" : "gcd";
        data.verbose		= verbose;
        data.async		= async;
	data.l1.x		= x;
	data.l1.y		= y;
        data.l2			= data.l1;	
        data.iter_count		= 0;
	data.log		= log;

	//copy clb array to start estimation with
	for(i=0;i<3;i++)
		starting_theta[i]=clb[i];

	//Optimizer initialization
        init_optimizer(&opt, clb, cub, 1e-5);

	//Create descriptors
        MORSE_Call(&data, *ncores,*gpus, *ts, *p_grid, *q_grid, *n,  0, 0);

        //Copy z to descriptor
	MORSE_Lapack_to_Tile( z, *n, data.descZ);

	//print summary message
        print_summary(1, *n, *ncores, *gpus, *ts, computation, 1, 1, 1);

	//main algorithm call
	START_TIMING(data.total_exec_time);
	nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
	nlopt_optimize(opt, starting_theta, &opt_f);
	STOP_TIMING(data.total_exec_time);

	//Print results to log files
	//print_result(&data, starting_theta, *n, 1, *ncores, *ts, 1, NULL, computation, 1, 1, data.final_loglik);

	//Destory descriptors & free memory
        nlopt_destroy(opt);
	MORSE_Desc_Destroy( &data.descC );
        MORSE_Desc_Destroy( &data.descZ );
        MORSE_Desc_Destroy( &data.descZcpy );
        MORSE_Desc_Destroy( &data.descproduct );
        MORSE_Desc_Destroy( &data.descdet );

        //copy local vector to global vector in R memory space
        for (j = 0; j < 3; j++)
                *(globalthetaout + j) = *(starting_theta + j);
}


void rexageostat_init(int *ncores, int *gpus, int *ts)
//! R-wrapper to initiate exageostat.
/*!  -- using dense or approximate computation
 * @param[in] ncores:           Pointer to the number of CPUs.
 * @param[in] gpus:             Pointer to the number of GPUs.
 * @param[in] ts:               Pointer to the tile size (MB) is used only in the case of HiCMA not MORSE.
 * */
{
	exageostat_init(ncores, gpus, ts);
}
void rexageostat_finalize()
//! R-wrapper to finalize exageostat.
{
        exageostat_finalize();
}

