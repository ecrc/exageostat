/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-14
 *
 **/
#include "../src/include/MLE.h"
#include "../r-wrappers/include/rwrappers.h"

int main(int argc, char **argv) {

 	//initialization
        char *theta;  //for testing case
        int n, ts, log, verbose;
        int i=0;
	char *dm;
	char *computation, *clb, *cub;
	int computation_int;
	int dm_int = 0;
        double time_opt = 0.0;
	int globalveclen;
        int gpus = 0, p_grid = 1, q_grid = 1, ncores = 1;
	int thetalen = 3;
        arguments arguments;
	double *vecs_out = NULL, *theta_out = NULL;
	double *initial_theta = NULL, *lb = NULL, *ub = NULL;
        double opt_tol = 1e-5;
        int opt_max_iters = 0;

        //Arguments default values
        set_args_default(&arguments);
        argp_parse (&argp, argc, argv, 0, 0, &arguments);

	//read inputs
	n		= atoi(arguments.N);
	ncores		= atoi(arguments.ncores);
        gpus		= atoi(arguments.gpus);
        p_grid		= atoi(arguments.p);
        q_grid		= atoi(arguments.q);    
	ts		= atoi(arguments.ts);    
	dm		= arguments.dm;
	computation	= arguments.computation; //approx or exact
	theta		= arguments.ikernel;
	clb		= arguments.olb;
	cub		= arguments.oub;
        opt_tol         = pow(10, -1.0 * atoi(arguments.opt_tol));
        opt_max_iters   = atoi(arguments.opt_max_iters);

        dm_int		= strcmp(dm, "ed") == 0 ? 0 : 1;
        computation_int = strcmp(computation, "exact") == 0 ? 0 : 1;
	globalveclen	= thetalen * n;
	//Memory allocation
	theta_out	= (double *) malloc(thetalen * sizeof(double));
	vecs_out        = (double *) malloc(globalveclen * sizeof(double));
	initial_theta   = (double *) malloc (thetalen * sizeof(double)); 
	lb           	= (double *) malloc(thetalen * sizeof(double));
        ub             	= (double *) malloc(thetalen * sizeof(double));
	//Parse inputs
	theta_parser2(initial_theta, theta);   
	theta_parser2(lb, clb);
        theta_parser2(ub, cub);


    rexageostat_init(&ncores,&gpus, &ts);
    rexageostat_gen_z(&n,   &ncores,   &gpus,  &ts,   &p_grid,  &q_grid,  &initial_theta[0],  &initial_theta[1],  &initial_theta[2],  &computation_int,  &dm_int, &globalveclen, vecs_out);
    rexageostat_likelihood(&n, &ncores, &gpus, &ts, &p_grid, &q_grid,  vecs_out, NULL,  &vecs_out[n], NULL,  &vecs_out[2*n], NULL,   lb, &thetalen, ub, thetalen, &computation_int, &dm_int, &opt_tol, &opt_max_iters, theta_out);
    rexageostat_finalize();
    
    fprintf(stderr,"Found Maximum at f(%g, %g, %g) \n", theta_out[0], theta_out[1], theta_out[2]);

    //free memory
    free(theta_out);
    free(vecs_out);
    free(ub);
    free(lb);
    return 0;
}

