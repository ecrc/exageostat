/**
 *
 * Copyright (c) 2017-2019  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE.c
 *
 * ExaGeoStat main functions.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2019-08-06
 *
 **/
#include "examples.h"
#include "../src/include/MLE.h"
#include "../r-wrappers/include/rwrappers.h"

int main(int argc, char **argv) {

 	//initialization
        char *theta;  //for testing case
        int n, dts, lts;// log, verbose;
        //int i=0;
	char *dm;
	char *computation, *clb, *cub;
	//int computation_int;
	int dm_int = 0, hicma_acc = 0, hicma_maxrank = 0, diag_thick = 0;
        //double time_opt = 0.0;
	int globalveclen;
        int gpus = 0, p_grid, q_grid, ncores;
	int thetalen = 3;
        char *kernel_fun;
        arguments arguments;
	double *vecs_out = NULL, *theta_out = NULL;
	double *initial_theta = NULL, *lb = NULL, *ub = NULL;
	double opt_tol;
	int opt_max_iters = 0;
	int num_params  = 0;
	char* locsFPath;
	location* locations;
	//Arguments default values
	set_args_default(&arguments);
	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	//read inputs
	//n		= atoi(arguments.N);   // n is determined based on the locations size.
	ncores		= atoi(arguments.ncores);
	gpus		= atoi(arguments.gpus);
	p_grid		= atoi(arguments.p);
	q_grid		= atoi(arguments.q);    
	diag_thick	= atoi(arguments.diag_thick);
	dts		= atoi(arguments.dts);
	lts             = atoi(arguments.lts); 
	hicma_acc	= atoi(arguments.acc);
	hicma_maxrank	= atoi(arguments.maxrank);   
	dm		= arguments.dm;
	computation	= arguments.computation; //approx or exact
	theta		= arguments.ikernel;
	clb		= arguments.olb;
	cub		= arguments.oub;
	opt_tol		= pow(10, -1.0 * atoi(arguments.opt_tol));
	opt_max_iters	= atoi(arguments.opt_max_iters);
	locsFPath       = arguments.locs_file;


	if(strcmp(kernel_fun, "univariate_matern_stationary")   == 0)
		num_params = 3;
	else if(strcmp(kernel_fun, "univariate_matern_nuggets_stationary")   == 0)
		num_params = 4;
	else if(strcmp(kernel_fun, "univariate_matern_non_stationary")   == 0)
		num_params = 9;
	else if(strcmp(kernel_fun, "bivariate_matern_flexible")   == 0)
		num_params = 13;
	else if(strcmp(kernel_fun, "bivariate_matern_parsimonious")   == 0)
		num_params = 6;
	else if(strcmp(kernel_fun, "bivariate_matern_parsimonious2")   == 0)
		num_params = 6;
	else if(strcmp(kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
		num_params = 6;
	else if(strcmp(kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
		num_params = 6;
        else if(strcmp(kernel_fun, "univariate_spacetime_matern_stationary")   == 0)
                num_params = 7;
	else
	{
		fprintf(stderr,"Choosen kernel is not exist!\n");
		fprintf("stderr, Called function is: %s\n",__func__);
		exit(0);
	}


	if(strcmp(locsFPath, "") == 0)
	{
		printf("Please determine the path of locations file using --locs_file argument\n");
		exit(0);
	}	


	//Read locations from a flat file.
	n = countlines(locsFPath);
	locations = readLocsFile(locsFPath, n);

	dm_int		= strcmp(dm, "ed") == 0 ? 0 : 1;

	globalveclen	=  n;

	//Memory allocation
	theta_out	= (double *) malloc((num_params+3) * sizeof(double));
	vecs_out        = (double *) malloc(globalveclen * sizeof(double));
	initial_theta   = (double *) malloc (thetalen * sizeof(double)); 
	lb           	= (double *) malloc(thetalen * sizeof(double));
	ub             	= (double *) malloc(thetalen * sizeof(double));

	//Parse inputs
	theta_parser2(initial_theta, theta, num_params);   
	theta_parser2(lb, clb, num_params);
	theta_parser2(ub, cub, num_params);



	rexageostat_init(&ncores,&gpus, &dts);
	//gen_z_exact(&n,   &ncores,   &gpus,  &dts,   &p_grid,  &q_grid,  &initial_theta[0],  &initial_theta[1],  &initial_theta[2], &dm_int, &globalveclen, vecs_out);

	gen_z_givenlocs_exact(&n, &ncores,  &gpus,
			&dts,  &p_grid, &q_grid,
			locations->x, NULL,  locations->y, NULL,
			&initial_theta[0],  &initial_theta[1],  &initial_theta[2],
			&dm_int, &globalveclen, vecs_out);


	if(strcmp (computation, "exact") == 0)
		mle_exact(&n, &ncores, &gpus,
				&dts, &p_grid, &q_grid,
				locations->x, NULL,  locations->y,
				NULL,  vecs_out, NULL,  
				lb, &thetalen, ub,
				&thetalen,  &dm_int, &opt_tol, 
				&opt_max_iters, theta_out);

	else if(strcmp (computation, "lr_approx") == 0)
		mle_tlr(&n, &ncores, &gpus,
				&lts, &p_grid, &q_grid, 
				locations->x, NULL,  locations->y,
				NULL,  vecs_out, NULL,   
				lb, &thetalen, ub, 
				&thetalen, &hicma_acc, &hicma_maxrank, 
				&dm_int, &opt_tol, &opt_max_iters, theta_out);	

	else if(strcmp (computation, "diag_approx") == 0)
		mle_dst(&n, &ncores, &gpus,
				&dts, &p_grid, &q_grid,
				locations->x, NULL,  locations->y, 
				NULL,  vecs_out, NULL, 
				lb, &thetalen, ub,
				&thetalen, &diag_thick,  &dm_int,
				&opt_tol, &opt_max_iters, theta_out);

	else if(strcmp (computation, "diag_approx") == 0)
		mle_dst(locations->x, NULL,  locations->y, NULL,  vecs_out, NULL, lb, &thetalen, ub, &thetalen, &diag_thick, &dm_int, &n, &opt_tol, &opt_max_iters, &ncores, &gpus, &lts, &p_grid, &q_grid,theta_out);
	// Old interfaces.
	//rexageostat_gen_z(&n,   &ncores,   &gpus,  &dts,   &p_grid,  &q_grid,  &initial_theta[0],  &initial_theta[1],  &initial_theta[2],  &computation_int,  &dm_int, &globalveclen, vecs_out);
	//rexageostat_likelihood(&n, &ncores, &gpus, &dts, &p_grid, &q_grid,  vecs_out, NULL,  &vecs_out[n], NULL,  &vecs_out[2*n], NULL,   lb, &thetalen, ub, &thetalen, &computation_int, &dm_int,  theta_out);

	rexageostat_finalize();
	printf("%f - %f - %f\n", theta_out[0], theta_out[1], theta_out[2]);

	//free memory
	free(theta_out);
	free(vecs_out);
	free(ub);
	free(lb);
	return 0;
}

