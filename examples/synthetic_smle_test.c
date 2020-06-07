/**
 *
 * Copyright (c) 2017-2019, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file synthetic_smle_test.c
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
	float *starting_theta;
	float *target_theta;
	float *initial_theta;  //for testing case
	int N, dts, lts, log;
	int i = 0;
	int seed = 0; 
	int zvecs = 1, nZmiss = 0, test = 0, gpus = 0, num_params = 3;
	int p_grid, q_grid, ncores;
	double  opt_f;
	arguments arguments;
	nlopt_opt opt;
	MLE_data data;
	location *locations = NULL;
	float prediction_error = 0.0;
	float* lb = (float *) malloc(3 * sizeof(float));
	float* up = (float *) malloc(3 * sizeof(float));
	int iseed[4]={seed, seed, seed, 1};

	//Arguments default values
	set_args_default(&arguments);
	opt=nlopt_create( NLOPT_LN_BOBYQA, 3);  //NLOPT_LN_BOBYQA  - NLOPT_GN_ORIG_DIRECT	
	argp_parse(&argp, argc, argv, 0, 0, &arguments);
	check_args(&arguments);

	//Memory allocation
	starting_theta	= (float *) malloc(3 * sizeof(float));
	initial_theta	= (float *) malloc(3 * sizeof(float));
	target_theta	= (float *) malloc(3 * sizeof(float));

	//MLE_data data initialization
	init(&test, &N,  &ncores, &gpus,
			&p_grid, &q_grid, &zvecs,
			&dts, &lts, &nZmiss,
			&log, initial_theta, starting_theta, 
			target_theta, lb, up, 
			&data, &arguments);

	//for single precision
	data.precision = 1;
	//printf("dts: %d, lts: %d\n", dts, lts);

	exageostat_init(&ncores, &gpus, &dts, &lts); 
	// Optimizater initialization
	//printf(" %d - %d \n", data.opt_tol, data.opt_max_iters);	
	init_optimizer(&opt, lb, up, pow(10, -1.0 * data.opt_tol));
	nlopt_set_maxeval(opt, data.opt_max_iters);

	//nomral random generation of e -- ei~N(0, 1) to generate Z
	double *Nrand	= (double *) malloc (N * zvecs * sizeof(double));
	float  *Nrands  = (float *) malloc (N * zvecs * sizeof(float));

	LAPACKE_dlarnv(3, iseed, N*zvecs, Nrand);
	//LAPACKE_slarnv(3, iseed, N*zvecs, Nrands);

	for(i=0; i < N*zvecs; i++)
		Nrands[i] = (float) Nrand[i];	
	free(Nrand);

	//data.l1.x=(float *) malloc(N * sizeof(float));
	//data.l1.y=(float *) malloc(N * sizeof(float));

	//uniform random generation for locations / read locations from disk
	//if(strcmp(data.locsFPath, "") == 0)
	locations = GenerateXYLoc(N, seed);
	//	else
	//		locations = readLocsFile(data.locsFPath, N);
	data.l1   = *locations;	

	int nZobs = strcmp(data.actualZFPath,"") == 0? (N-nZmiss) : N;

	if(strcmp (data.computation, "exact") == 0)
	{
		if(data.precision == 0)
			MORSE_dmle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N,  nZobs, nZmiss);
		else
			MORSE_smle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N,  nZobs, nZmiss);		
	}
	else if (strcmp (data.computation, "diag_approx") == 0)
		MORSE_dmle_diag_Call(&data, ncores, gpus, dts, p_grid, q_grid, N,  nZobs, nZmiss);		
#if defined(EXAGEOSTAT_USE_HICMA)
	else if (strcmp (data.computation, "lr_approx") == 0)
	{
		HICMA_dmle_Call(&data, ncores, gpus, lts, p_grid, q_grid, N,  nZobs, nZmiss);
		if (test == 1)
		{
			data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT;
			printf("test case: HICMA_STARSH_PROB_GEOSTAT\n");
		}
		else
			data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
	}
#endif		


	for(i=0;i<zvecs;i++)
	{

		print_summary(test, N, ncores, gpus, dts, lts,  data.computation, zvecs, p_grid, q_grid, data.precision);

		if(arguments.profile == 1)
		{
			starpu_fxt_autostart_profiling(0);
			starpu_fxt_start_profiling();
		}

		//MLE_zvg(&data, &Nrand[i*N], initial_theta, N, dts, log) ;
		MLE_zvg(&data, &Nrands[i*N], initial_theta, N, dts, log, p_grid, q_grid);
		if(log == 1 && test == 1)
			init_log(&data);                 

		START_TIMING(data.total_exec_time);
		nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
		nlopt_optimize(opt, starting_theta, &opt_f);
		STOP_TIMING(data.total_exec_time);




		if(nZmiss != 0){

			//initialization
			float *Zobs;
			float *Zactual;
			float *Zmiss;
			int i=0;

			for (i=0;i<10;i++)

				//memory allocation
				Zobs 	= (float *) malloc(nZobs * sizeof(float));
			Zactual	= (float *) malloc(nZmiss * sizeof(float));
			Zmiss	= (float *) malloc(nZmiss * sizeof(float));
			if(strcmp (data.computation, "exact") == 0 || strcmp (data.computation, "diag_approx") == 0)
				prediction_init(&data, nZmiss, nZobs, dts, p_grid, q_grid, 1);	
#if defined(EXAGEOSTAT_USE_HICMA)
			else if (strcmp (data.computation, "lr_approx") == 0)
				prediction_init(&data, nZmiss, nZobs, lts, p_grid, q_grid, 1);						
#endif
			int j = 0;
			for (j = 0; j < 1; j++)
			{
				//printf("j = %d\n",j);
				pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);
				//generate_interior_points(&data, Zobs, NULL, nZmiss, nZobs, N);
				if (strcmp (data.computation, "exact") == 0)
					prediction_error = MORSE_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N);
				else if (strcmp (data.computation, "diag_approx") == 0)
					prediction_error = MORSE_dmle_diag_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N);
#if defined(EXAGEOSTAT_USE_HICMA)
				else if (strcmp (data.computation, "lr_approx") == 0)
				{
					data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
					prediction_error = HICMA_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual, Zmiss, N, lts);
				}
#endif
				int index=0;
				for (index=0; index< nZmiss; index++)
					printf ("(%f, %f)\n ", Zactual[index], Zmiss[index]);

				fprintf(stderr,"Prediction Error: %f \n", prediction_error);
				write_to_thetafile("theta-pred.txt", starting_theta, num_params, prediction_error, j, 0, 0, 0 ,0);

			}



			prediction_finalize(&data);
			//free memory
			free(Zactual);
			free(Zobs);
		}
		print_result(&data, starting_theta, N, zvecs, ncores, dts, test, arguments.ikernel, data.computation, p_grid, q_grid, data.final_loglik, prediction_error);
		if(log == 1 && test == 1)
			finalize_log(&data);

	}

	nlopt_destroy(opt);
	MLE_Finalize(&data);

	if(arguments.profile == 1)
	{
		starpu_fxt_stop_profiling();
		RUNTIME_profiling_display_efficiency();
	}

	return 0;
}


