/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include "../src/include/MLE.h"

int main(int argc, char **argv) {

        //initialization
        double * starting_theta;
        double * target_theta;
        double * initial_theta;  //for testing case
        int N, ts, log;
        int i = 0;
        int zvecs = 1, nZmiss = 0, test = 0, gpus = 0;
        int p_grid, q_grid, ncores;
        double  opt_f;
        arguments arguments;
        nlopt_opt opt;
        double * Nrand;
        MLE_data data;
        double* lb = (double *) malloc(3 * sizeof(double));
	double* up = (double *) malloc(3 * sizeof(double));
	int iseed[4]={0, 0, 0, 1};

        //Arguments default values
        set_args_default(&arguments);
        argp_parse(&argp, argc, argv, 0, 0, &arguments);
	check_args(&arguments);

	//Memory allocation
        starting_theta	= (double *) malloc(3 * sizeof(double));
        initial_theta	= (double *) malloc(3 * sizeof(double));
        target_theta	= (double *) malloc(3 * sizeof(double));

        //MLE_data data initialization
        init(&test, &N,  &ncores, &gpus, &p_grid, &q_grid, &zvecs, &ts, &nZmiss, &log, initial_theta, starting_theta, target_theta, lb, up, &data, &arguments);
       
	exageostat_init(&ncores, &gpus, &ts); 
	
	//kernel parsing
        init_optimizer(&opt, lb, up, 1e-5);


        //nomral random generation of e -- ei~N(0, 1) to generate Z
        Nrand	= (double *) malloc (N * zvecs * sizeof(double));
        LAPACKE_dlarnv(3, iseed, N*zvecs, Nrand);

	data.l1.x=(double *) malloc(N * sizeof(double));
        data.l1.y=(double *) malloc(N * sizeof(double));

        if(strcmp (data.computation, "exact") == 0)	
        	//uniform random generation for locations / read locations from disk
		GenerateXYLoc(N, data.locsFPath, &data.l1);


	int nZobs = strcmp(data.actualZFPath,"") == 0? (N-nZmiss) : N;

	if(strcmp (data.computation, "exact") == 0)
		MORSE_Call(&data, ncores, gpus, ts,p_grid, q_grid, N,  nZobs, nZmiss);


        for(i=0;i<zvecs;i++)
        {

                print_summary(test, N, ncores, gpus, ts, data.computation, zvecs, p_grid, q_grid);

		if(arguments.profile == 1)
		{
			starpu_fxt_autostart_profiling(0);
			starpu_fxt_start_profiling();
		}
		
                MLE_zvg(&data, &Nrand[i*N], initial_theta, N,ts, test, log) ;


	
		if(log == 1 && test == 1)
			init_log(&data);                 

                START_TIMING(data.total_exec_time);
                nlopt_set_max_objective(opt, MLE_alg, (void *)&data);
                nlopt_optimize(opt, starting_theta, &opt_f);
                STOP_TIMING(data.total_exec_time);
              

		print_result(&data, starting_theta, N, zvecs, ncores, ts, test, arguments.ikernel, data.computation, p_grid, q_grid, data.final_loglik);


	        if(nZmiss != 0){
			double prediction_error = MORSE_MLE_Predict_Tile(&data, starting_theta, nZmiss, nZobs, N);
	                fprintf(stderr,"Prediction Error: %f \n",prediction_error);
        	}

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


