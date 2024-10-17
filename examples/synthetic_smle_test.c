/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
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
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "examples.h"
#include "../src/include/MLE.h"

int main(int argc, char** argv) {

    //initialization
    double* starting_theta;
    double* target_theta;
    double* initial_theta;  //for testing case
    int num_params = 0;
    int N, dts, lts, log;
    int i = 0;
    int j = 0;
    int p = 0;         //univariate/bivariate/multivariate
    int zvecs = 1, nZmiss = 0, test = 0, gpus = 0;
    double x_max, x_min, y_max, y_min;
    int p_grid, q_grid, ncores;
    double opt_f;
    arguments arguments;
    double all_time = 0.0;
    double pred_time = 0.0;
    nlopt_opt opt;
    MLE_data data;
    int seed = 0;
    location *locations = NULL;
    double prediction_error = 0.0;
    //Arguments default values
    set_args_default(&arguments);
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    check_args(&arguments);

    if (strcmp(arguments.kernel_fun, "univariate_matern_stationary") == 0) {
        num_params = 3;
        p = 1;
    } else if (strcmp(arguments.kernel_fun, "univariate_matern_nuggets_stationary") == 0) {
        num_params = 4;
        p = 1;
    } else if (strcmp(arguments.kernel_fun, "univariate_matern_non_stationary") == 0) {
        num_params = 9;
        p = 1;
    } else if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
               strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0) {
        num_params = 6;
        p = 1;
    } else if (strcmp(arguments.kernel_fun, "bivariate_matern_flexible") == 0) {
        num_params = 11;
        p = 2;
    } else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious") == 0 ||
               strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious_profile") == 0) {
        num_params = 6;
        p = 2;
    } else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
               strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        num_params = 6;
        p = 2;
    } else if (strcmp(arguments.kernel_fun, "univariate_spacetime_matern_stationary") == 0) {
        num_params = 7;
        p = 1;
    } else if (strcmp(arguments.kernel_fun, "bivariate_spacetime_matern_stationary") == 0) {
        num_params = 10;
        p = 2;
    } else if (strcmp(arguments.kernel_fun, "trivariate_matern_parsimonious") == 0 ||
               strcmp(arguments.kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {
        num_params = 10;
        p = 3;
    } else if (strcmp(arguments.kernel_fun, "univariate_matern_non_stat") == 0) {
        num_params = 8;
        p = 1;
    } else {
        fprintf(stderr, "Choosen kernel is not exist(20)!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }


    double* lb = (double* ) malloc(num_params * sizeof(double));
    double* up = (double* ) malloc(num_params * sizeof(double));
    int iseed[4] = {seed, seed, seed, 1};

    //Memory allocation
    starting_theta = (double* ) malloc(num_params * sizeof(double));
    initial_theta = (double* ) malloc(num_params * sizeof(double));
    target_theta = (double* ) malloc(num_params * sizeof(double));

    //MLE_data data initialization
    init(&test, &N, &ncores,
         &gpus, &p_grid, &q_grid,
         &zvecs, &dts, &lts,
         &nZmiss, &log, initial_theta,
         starting_theta, target_theta, lb,
         up, &data, &arguments);

    if (strcmp(arguments.kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        p = data.time_slots;
    if (strcmp(arguments.kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        p = 2 * data.time_slots;
    if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0
        || strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2_profile") == 0
            ) {
        if (N % dts != 0) {
            printf("please use N divisible by dts, only with parsimonious2\n");
            exit(0);
        }

    }

    int nZobs = strcmp(data.actualZFPath, "") == 0 ? (N - nZmiss) : N;
    //To support multivariate case
    N = p * N;
    results.problem_size = N;
    exageostat_init(&ncores, &gpus, &dts, &lts);
    // Optimizater initialization
    //NLOPT_LN_BOBYQA && NLOPT_LN_SBPLX
    printf("num_params: %d,\n", num_params);
    opt = nlopt_create(NLOPT_LN_BOBYQA, num_params);
    init_optimizer(&opt, lb, up, pow(10, -1.0 * data.opt_tol));
    nlopt_set_maxeval(opt, data.opt_max_iters);

    //for single precision
    data.precision = 1;
    //nomral random generation of e -- ei~N(0, 1) to generate Z
    double* Nrand_d = (double* ) malloc(N * zvecs * sizeof(double));
    float* Nrand = (float* ) malloc(N * zvecs * sizeof(float));

    LAPACKE_dlarnv(3, iseed, N * zvecs, Nrand_d);
    //LAPACKE_slarnv(3, iseed, N*zvecs, Nrands);

    for (i = 0; i < N * zvecs; i++)
        Nrand[i] = (float) Nrand_d[i];
    free(Nrand_d);

    printf("%s\n", arguments.dim);

    //uniform random generation for locations / read locations from disk
    if (strcmp(arguments.dim, "2d") == 0) {
        printf("2D example......\n");
        locations = GenerateXYLoc(N / p, seed);
        data.l1.z = NULL;
    } else if (strcmp(arguments.dim, "3d") == 0) {
        printf("3D example......\n");
        //// TODO: FIX THIS TOO.
        locations = GenerateXYZLoc(N / p, seed);
    } else if (strcmp(arguments.dim, "st") == 0) {
        printf("ST example......\n");
        locations = GenerateXYLoc_ST(N / p, data.time_slots, seed);
    } else {
        printf("Input dimension is not supported. Please use 2d, 3d, or st\n\n");
        exit(0);
    }
    data.l1 = *locations;
    //find x_max, x_min, y_max, y_min
    if (strcmp(data.kernel_fun, "univariate_matern_non_stationary") == 0) {
        x_max = data.l1.x[0];
        x_min = data.l1.x[0];
        y_max = data.l1.y[0];
        y_min = data.l1.y[0];
        for (i = 1; i < N; i++) {
            if (x_max < data.l1.x[i])
                x_max = data.l1.x[i];
            if (x_min > data.l1.x[i])
                x_min = data.l1.x[i];
            if (y_max < data.l1.y[i])
                y_max = data.l1.y[i];
            if (y_max > data.l1.y[i])
                y_max = data.l1.y[i];
        }

        data.lm.x[0] = x_min + (x_max - x_min) / 2;
        data.lm.y[0] = y_min + (y_max - y_min) / 2;
        printf(" The central point is ( %f, %f)\n", data.lm.x[0], data.lm.y[0]);
    }
    //***************************


    if (strcmp(data.computation, "exact") == 0) {
        if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
            strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0) {
            EXAGEOSTAT_dmle_ng_Call(&data, ncores, gpus, dts, p_grid, q_grid, N, nZobs, nZmiss);
        } else {
            if (data.precision == 0)
                EXAGEOSTAT_dmle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N, nZobs, nZmiss);
            else
                EXAGEOSTAT_smle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N, nZobs, nZmiss);
        }
    } else if (strcmp(data.computation, "diag_approx") == 0)
        EXAGEOSTAT_dmle_diag_Call(&data, ncores, gpus, dts, p_grid, q_grid, N, nZobs, nZmiss);
#if defined(EXAGEOSTAT_USE_HICMA)
    else if (strcmp(data.computation, "lr_approx") == 0) {
        if (test == 1) {
            if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0)
                EXAGEOSTAT_TLR_dmle_ng_Call(&data, ncores, gpus, lts, p_grid, q_grid, N, nZobs, nZmiss);
            else
                EXAGEOSTAT_TLR_dmle_Call(&data, ncores, gpus, lts, dts, p_grid, q_grid, N, nZobs, nZmiss);

            if (strcmp(arguments.kernel_fun, "univariate_matern_stationary") == 0 ||
                strcmp(arguments.kernel_fun, "univariate_matern_nuggets_stationary") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT;
            else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE;
            else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE;
            if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN;
        } else {

            if (strcmp(arguments.kernel_fun, "univariate_matern_stationary") == 0 ||
                strcmp(arguments.kernel_fun, "univariate_matern_nuggets_stationary") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
            else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT;
            else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT;
            if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0)
                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN_POINT;
        }
    }
#endif


    for (i = 0; i < zvecs; i++) {
        //output
        results.z_sample = i;
        print_summary(test, N, ncores, gpus, dts, lts, data.computation, zvecs, p_grid, q_grid, data.precision);

        //reset starting_theta vector
        for (j = 0; j < num_params; j++)
            starting_theta[j] = lb[j];

        if (arguments.profile == 1) {
            starpu_fxt_autostart_profiling(0);
            starpu_fxt_start_profiling();
        }

        //MLE_zvg(&data, &Nrand[i*N], initial_theta, N, dts, log) ;
        if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
            strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0) {
            MLE_ng_zvg(&data, &Nrand[i * N], initial_theta, N, dts, log, p_grid, q_grid);
            if (log == 1 && test == 1)
                init_log(&data);
            //use a starting point for g and h assuming Gaussian data
            starting_theta[4] = 0.2;
            starting_theta[5] = 0.2;

            START_TIMING(data.total_exec_time);
            nlopt_set_max_objective(opt, MLE_ng_alg, (void *) &data);
            nlopt_optimize(opt, starting_theta, &opt_f);
            STOP_TIMING(data.total_exec_time);
            for (j = 0; j < num_params; j++) {
                results.estimated_theta[j] = starting_theta[j];
            }
            results.final_loglik = opt_f;
        } else {

            MLE_zvg(&data, &Nrand[i * N], initial_theta, N, dts, log, p_grid, q_grid);
            if (log == 1 && test == 1)
                init_log(&data);

            START_TIMING(data.total_exec_time);
            nlopt_set_max_objective(opt, MLE_alg, (void *) &data);
            nlopt_optimize(opt, starting_theta, &opt_f);
            STOP_TIMING(data.total_exec_time);
        }

        //Calculate the Fisher matrix
        if (data.fisher == 1) {
            //print the fisher matrix
            int k = 0;
            for (k = 0; k < 9; k++) {
                if (k % 3 == 0)
                    printf("\n");
            }
            fprintf(stderr, "\n");
        }
        if (nZmiss != 0) {

            //initialization
            double* Zobs;
            double* Zactual;
            double* Zmiss;

            int i = 0;
            double avg_pred_value = 0.0;
            double avg_pred_value1 = 0.0;
            double avg_pred_value2 = 0.0;
            double avg_pred_value3 = 0.0;
            int pred_samples = 1;

            if (data.mloe_mmom == 1 || data.mloe_mmom_async == 1) {
                nZobs = strcmp(data.actualZFPath, "") == 0 ? (N / p - nZmiss) : N;
#if defined(EXAGEOSTAT_USE_HICMA)
                if (strcmp(data.computation, "lr_approx") == 0) {
                    if (strcmp(arguments.kernel_fun, "univariate_matern_stationary") == 0 ||
                        strcmp(arguments.kernel_fun, "univariate_matern_nuggets_stationary") == 0)
                        data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
                    else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious") == 0)
                        data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT;
                    else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0)
                        data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT;
                    if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                        strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0)
                        data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN_POINT;
                }
#endif
                printf("NZobs=%d, nZmiss=%d, N=%d,\n", nZobs, nZmiss, N);
                //exit(0);
                Zobs = (double* ) malloc(nZobs * sizeof(double));
                Zactual = (double* ) malloc(nZmiss * sizeof(double));
                if (p == 2)
                    pick_random_points_noshuffle(&data, Zobs, Zactual, nZmiss, nZobs, N);

                else
                    pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);
                mloe_mmom_init(&data, nZmiss, nZobs, dts, p_grid, q_grid);
                START_TIMING(all_time);
                if (data.mloe_mmom == 1)
                    EXAGEOSTAT_dmle_mloe_mmom_Tile(&data, initial_theta, starting_theta, nZmiss, nZobs, N);
                if (data.mloe_mmom_async == 1)
                    EXAGEOSTAT_dmle_mloe_mmom_Tile_Async(&data, initial_theta, starting_theta, nZmiss, nZobs, N);
                //TO BE REMOVED
                data.kernel_fun = arguments.kernel_fun;
                STOP_TIMING(all_time);
                free(Zobs);
                free(Zactual);
                fprintf(stderr, " ---- mloe_mmom Time(main): %6.2f seconds\n\n", all_time);

                MLOE_MMOM_Finalize(&data);
            }

            if (data.idw == 1) {
                //memory allocation

                Zobs = (double* ) malloc(p * nZobs * sizeof(double));
                Zactual = (double* ) malloc(p * nZmiss * sizeof(double));
                Zmiss = (double* ) malloc(p * nZmiss * sizeof(double));
                double* mspe = (double* ) malloc(4 * sizeof(double));
                int j = 0;
                for (j = 0; j < pred_samples; j++) {
                    printf("nZobs = %d\n", p * nZobs);
                    if (p == 2)
                        pick_random_points2(&data, Zobs, Zactual, nZmiss, nZobs, N);
                    else if (p == 3)
                        pick_random_points3(&data, Zobs, Zactual, nZmiss, nZobs, N);
                    else
                        pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);

                    START_TIMING(pred_time);
                    mspe = pred_idw(&data, Zmiss, Zactual, Zobs, nZmiss, nZobs);
                    STOP_TIMING(pred_time);
#if defined(CHAMELEON_USE_MPI)
                    if(CHAMELEON_My_Mpi_Rank() == 0)
                    {
#endif
                    int index = 0;
                    for (index = 0; index < nZmiss; index++)
                        printf("(%3.6f, %3.6f)\n ", Zactual[index], Zmiss[index]);
                    fprintf(stderr, "Prediction Error (IDW): %3.9f  -  %3.9f -  %3.9f\n", mspe[0], mspe[1], mspe[2]);
#if defined(CHAMELEON_USE_MPI)
                    }
#endif
                    avg_pred_value += mspe[0];
                    avg_pred_value1 += mspe[1];
                    avg_pred_value2 += mspe[2];
                    avg_pred_value3 += mspe[3];
                }

#if defined(CHAMELEON_USE_MPI)
                if(CHAMELEON_My_Mpi_Rank() == 0)
                {
#endif
                fprintf(stderr, "Average prediction Error (IDW): %3.9f  -  %3.9f -  %3.9f\n", avg_pred_value,
                        avg_pred_value1, avg_pred_value2);

#if defined(CHAMELEON_USE_MPI)
                }
#endif

                //free memory
                free(Zactual);
                free(Zobs);
                free(Zmiss);
            }
            if (data.mspe == 1) {
                //memory allocation
                Zobs = (double* ) malloc(p * nZobs * sizeof(double));
                Zactual = (double* ) malloc(p * nZmiss * sizeof(double));
                Zmiss = (double* ) malloc(p * nZmiss * sizeof(double));
                if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                    strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0) {
                    if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "diag_approx") == 0)
                        EXAGEOSTAT_dmle_ng_Predict_Allocate(&data, nZmiss, nZobs, dts, p_grid, q_grid, 1);
#if defined(EXAGEOSTAT_USE_HICMA)
                    else if (strcmp(data.computation, "lr_approx") == 0)
                        EXAGEOSTAT_TLR_dmle_ng_Predict_Allocate(&data, nZmiss, nZobs, lts, p_grid, q_grid, 1);
#endif
                } else {
                    if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "diag_approx") == 0)
                        prediction_init(&data, nZmiss, nZobs, dts, p_grid, q_grid, 1);
#if defined(EXAGEOSTAT_USE_HICMA)
                    else if (strcmp(data.computation, "lr_approx") == 0)
                        prediction_init(&data, nZmiss, nZobs, lts, p_grid, q_grid, 1);
#endif
                }

                int j = 0;
                for (j = 0; j < pred_samples; j++) {
                    printf("nZobs = %d\n", p * nZobs);
                    if (p == 2)
                        pick_random_points2(&data, Zobs, Zactual, nZmiss, nZobs, N);
                    else if (p == 3)
                        pick_random_points3(&data, Zobs, Zactual, nZmiss, nZobs, N);
                    else
                        pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);

                    START_TIMING(pred_time);
                    if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                        strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0) {
                        if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "diag_approx") == 0)
                            prediction_error = EXAGEOSTAT_dmle_ng_Predict_Tile(&data, starting_theta, p * nZmiss,
                                                                              p * nZobs, Zobs, Zactual, Zmiss, N);
#if defined(EXAGEOSTAT_USE_HICMA)
                        else if (strcmp(data.computation, "lr_approx") == 0) {
                            if (strcmp(arguments.kernel_fun, "univariate_matern_stationary") == 0 ||
                                strcmp(arguments.kernel_fun, "univariate_matern_nuggets_stationary") == 0)
                                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
                            else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious") == 0)
                                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT;
                            else if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0)
                                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT;
                            if (strcmp(arguments.kernel_fun, "univariate_matern_non_gaussian") == 0 ||
                                strcmp(arguments.kernel_fun, "univariate_exp_non_gaussian") == 0)
                                data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN_POINT;


                            prediction_error = EXAGEOSTAT_TLR_dmle_ng_Predict_Tile(&data, starting_theta, p * nZmiss, p * nZobs,
                                                                          Zobs, Zactual, Zmiss, N, lts);
                        }
#endif

                    } else {
                        if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "lr_approx") == 0)
                            prediction_error = EXAGEOSTAT_dmle_Predict_Tile(&data, starting_theta, p * nZmiss, p * nZobs,
                                                                           Zobs, Zactual, Zmiss, N);
                        else if (strcmp(data.computation, "diag_approx") == 0)
                            prediction_error = EXAGEOSTAT_dmle_diag_Predict_Tile(&data, starting_theta, nZmiss, nZobs,
                                                                                Zobs, Zactual, Zmiss, N);
                    }

                    STOP_TIMING(pred_time);
#if defined(CHAMELEON_USE_MPI)
                    if(CHAMELEON_My_Mpi_Rank() == 0)
                    {
#endif
                    int index = 0;
                    for (index = 0; index < nZmiss; index++)
                        printf("(%3.6f, %3.6f)\n ", Zactual[index], Zmiss[index]);

                    fprintf(stderr, "Prediction Error: %3.9f \n", prediction_error);
#if defined(CHAMELEON_USE_MPI)
                    }
#endif
                    avg_pred_value += prediction_error;
                    avg_pred_value1 += data.mserror1;
                    avg_pred_value2 += data.mserror2;
                    avg_pred_value3 += data.mserror3;
                }


                prediction_finalize(&data);
                //free memory
                free(Zactual);
                free(Zobs);
                free(Zmiss);
            }
            char buf[30];
            char str[80];
            strcpy(str, arguments.kernel_fun);
            for (int i = 0; i < num_params; i++) {

                sprintf(buf, "%0.3f-", initial_theta[i]);
                strcat(str, buf);
            }
            strcat(str, data.computation);
            sprintf(buf, "%0.0f-", data.hicma_acc);
            strcat(str, buf);
            strcat(str, "-theta.txt");
            if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious_profile") == 0) {
                starting_theta[0] = data.variance1;
                starting_theta[1] = data.variance2;
            } else if (strcmp(arguments.kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {
                starting_theta[0] = data.variance1;
                starting_theta[1] = data.variance2;
                starting_theta[2] = data.variance3;
            }

            write_to_estimatedtheta(str, starting_theta, num_params, N / p, pred_time, all_time,
                                    (avg_pred_value1 /= pred_samples), (avg_pred_value2 /= pred_samples),
                                    (avg_pred_value3 /= pred_samples), (avg_pred_value /= pred_samples), data.mloe,
                                    data.mmom, zvecs);
        }

        print_result(&data, starting_theta, N, zvecs, ncores, lts, test, initial_theta, data.computation, p_grid,
                     q_grid, data.final_loglik, prediction_error);

        if (log == 1 && test == 1)
            finalize_log(&data);

    }

    nlopt_destroy(opt);

    if (arguments.profile == 1) {
        starpu_fxt_stop_profiling();
        RUNTIME_profiling_display_efficiency();
    }


    CHAMELEON_Finalize();
    return 0;
}