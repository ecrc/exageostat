/**
 *
 * Copyright (c) 2017-2019  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file zgen_mle_test.c
 *
 * A complete example to test ExaGeoStat supported function (i.e., dataset generator, Maximum Likelihood Function (MLE), Prediction)
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2019-08-06
 *
 **/
#include "examples.h"
#include "../src/include/MLE.h"

int main(int argc, char **argv) {

        //initialization
        double * starting_theta;
        double * target_theta;
        double * initial_theta;  //for testing case
        int N, lts, dts, log;
        //int i = 0;
        int zvecs = 1, nZmiss = 0, test = 0, gpus = 0;
        int p_grid, q_grid, ncores;
        double  opt_f;
        arguments arguments;
        nlopt_opt opt;
        double *streamdata;
        MLE_data data;
	location *locations;
	location *missing_locations;
        double* lb = (double *) malloc(3 * sizeof(double));
	double* up = (double *) malloc(3 * sizeof(double));
	//int iseed[4]={0, 0, 0, 1};

        //Arguments default values
        set_args_default(&arguments);
        argp_parse(&argp, argc, argv, 0, 0, &arguments);
	check_args(&arguments);

	//Memory allocation
        starting_theta	= (double *) malloc(3 * sizeof(double));
        initial_theta	= (double *) malloc(3 * sizeof(double));
        target_theta	= (double *) malloc(3 * sizeof(double));

        //MLE_data data initialization
        init(&test, &N,  &ncores, &gpus, &p_grid, &q_grid, &zvecs, &dts, &lts, &nZmiss, &log, initial_theta, starting_theta, target_theta, lb, up, &data, &arguments);
       
	exageostat_init(&ncores, &gpus, &dts, &lts); 
	
	//kernel parsing
        opt=nlopt_create( NLOPT_LN_BOBYQA, 3);  //NLOPT_LN_BOBYQA  - NLOPT_GN_ORIG_DIRECT
        init_optimizer(&opt, lb, up, pow(10, -1.0 * data.opt_tol));
        nlopt_set_maxeval(opt, data.opt_max_iters);


	//data.l1.x=(double *) malloc(N * sizeof(double));
        //data.l1.y=(double *) malloc(N * sizeof(double));

       	//Read locations from a flat file.
	N = countlines(data.locsFPath);
	locations = readLocsFile(data.locsFPath, N);
	data.l1   = *locations;	

	int nZobs = strcmp(data.actualZFPath,"") == 0? (N-nZmiss) : N;

	if(strcmp (data.computation, "exact") == 0)
		MORSE_dmle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N,  nZobs, nZmiss);
	else if (strcmp (data.computation, "diag_approx") == 0)
                MORSE_dmle_diag_Call(&data, ncores, gpus, dts, p_grid, q_grid, N,  nZobs, nZmiss);		
        #if defined(EXAGEOSTAT_USE_HICMA)
	else if (strcmp (data.computation, "lr_approx") == 0)
                {
                        HICMA_dmle_Call(&data, ncores, gpus, lts, p_grid, q_grid, N,  nZobs, nZmiss);
                        if (test == 1)
                                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT;
                        else
                                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
                }
	#endif

	print_summary(test, N, ncores, gpus, dts, lts,  data.computation, zvecs, p_grid, q_grid, data.precision);

	if(arguments.profile == 1)
	{
		starpu_fxt_autostart_profiling(0);
		starpu_fxt_start_profiling();
	}
	//read observation file	
	streamdata = readObsFile(data.obsFPath, N);
	locations_obs_zsort_inplace(N, locations, streamdata);
	if(strcmp (data.computation, "exact") == 0 || strcmp (data.computation, "diag_approx") == 0)
		MORSE_MLE_dzcpy(&data, streamdata);
#if defined(EXAGEOSTAT_USE_HICMA)
	else if (strcmp (data.computation, "lr_approx") == 0)	
		HICMA_MLE_zcpy(&data, streamdata);				
#endif
	if(log == 1 && test == 1)
		init_log(&data);                 

	START_TIMING(data.total_exec_time);
	nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
	nlopt_optimize(opt, starting_theta, &opt_f);
	STOP_TIMING(data.total_exec_time);


	print_result(&data, starting_theta, N, zvecs, ncores, dts, test, arguments.ikernel, data.computation, p_grid, q_grid, data.final_loglik);

	if (strcmp(data.actualZLocFPath,"") != 0)
	{
		printf( "%s ========\n", data.actualZLocFPath);
		nZmiss = countlines(data.actualZLocFPath);
		missing_locations = readLocsFile(data.actualZLocFPath, N);
	}	

	if(nZmiss != 0){

		//initialization
		double *Zobs;
		double *Zactual;
		double *Zmiss;

		//memory allocation
		Zobs    = (double *) malloc(nZobs * sizeof(double));
		Zactual = (double *) malloc(nZmiss * sizeof(double));
		Zmiss   = (double *) malloc(nZmiss * sizeof(double));

		if(strcmp (data.computation, "exact") == 0 || strcmp (data.computation, "diag_approx") == 0)
			prediction_init(&data, nZmiss, nZobs, dts, p_grid, q_grid, 1);
#if defined(EXAGEOSTAT_USE_HICMA)
		else if (strcmp (data.computation, "lr_approx") == 0)
			prediction_init(&data, nZmiss, nZobs, lts, p_grid, q_grid, 1);
#endif

		int j = 0;
		for (j = 0; j < 1; j++)
		{
			if (strcmp(data.actualZLocFPath,"") == 0)
				pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);
			else
			{
				Zactual = readObsFile(data.actualZFPath, nZmiss);
				MLE_get_zobs(&data, Zobs, N);
				data.lmiss = *missing_locations;
				data.lobs  = *locations;
			}
			// generate_interior_points(&data, Zobs, NULL, nZmiss, nZobs, N);
			double prediction_error = 0.0;
			if (strcmp (data.computation, "exact") == 0)
				prediction_error = MORSE_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N);
			else if (strcmp (data.computation, "diag_approx") == 0)
				prediction_error = MORSE_dmle_diag_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N);
#if defined(EXAGEOSTAT_USE_HICMA)
			else if (strcmp (data.computation, "lr_approx") == 0)
				prediction_error = HICMA_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N, lts);
#endif
			fprintf(stderr,"Prediction Error: %f \n", prediction_error);
		}

		int index=0;
		for (index=0; index< nZmiss; index++)
			printf ("(%f, %f)\n ", Zactual[index], Zmiss[index]);

		prediction_finalize(&data);
		//free memory
		free(Zactual);
		free(Zobs);
		free(Zmiss);
	}

	if(log == 1 && test == 1)
		finalize_log(&data);





	nlopt_destroy(opt);
	MLE_Finalize(&data);

	if(arguments.profile == 1)
	{
		starpu_fxt_stop_profiling();
		RUNTIME_profiling_display_efficiency();
	}

	return 0;
}


