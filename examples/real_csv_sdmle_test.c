/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file real_csv_sdmle_test.c
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
    int N, lts, dts, log;
    int zvecs = 1, nZmiss = 0, test = 0, gpus = 0, num_params = 3;
    int p_grid, q_grid, ncores;
    double opt_f;
    int p = 1;
    double all_time = 0.0;
    double pred_time = 0.0;
    arguments arguments;
    nlopt_opt opt;
    double* streamdata;
    double* streamdata1;
    double* streamdata2;
    double* streamdata3;
    MLE_data data;
    location *locations;
    location *missing_locations;
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
    }
    double* lb = (double* ) malloc(num_params * sizeof(double));
    double* up = (double* ) malloc(num_params * sizeof(double));


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

    if (strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2") == 0
        || strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious_profile") == 0
        || strcmp(arguments.kernel_fun, "bivariate_matern_parsimonious2_profile") == 0
            ) {
        if (N % dts != 0) {
            printf("please use N divisible by dts, only with parsimonious2\n");
            exit(0);
        }

    }

    //kernel parsing
    opt = nlopt_create(NLOPT_LN_BOBYQA, num_params);  //NLOPT_LN_BOBYQA  - NLOPT_GN_ORIG_DIRECT
    init_optimizer(&opt, lb, up, pow(10, -1.0 * data.opt_tol));
    nlopt_set_maxeval(opt, data.opt_max_iters);


    data.precision = 2;

    //Read locations from a flat file.
    N = countlines(data.locsFPath);
    locations = readLocsFile(data.locsFPath, N);
    data.l1 = *locations;

    int nZobs = strcmp(data.actualZFPath, "") == 0 ? (N - nZmiss) : N;
    //To support multivariate case
    N = p * N;

    exageostat_init(&ncores, &gpus, &dts, &lts);

    if (strcmp(data.computation, "exact") == 0) {
        EXAGEOSTAT_sdmle_Call(&data, ncores, gpus, dts, p_grid, q_grid, N, nZobs, nZmiss);
        EXAGEOSTAT_MLE_sdregister_Tile(&data);

    }

    print_summary(test, N, ncores, gpus, dts, lts, data.computation, zvecs, p_grid, q_grid, data.precision);

    if (arguments.profile == 1) {
        starpu_fxt_autostart_profiling(0);
        starpu_fxt_start_profiling();
    }
    //read observation file
    if (p == 1) {
        //for wind dataset
        streamdata = (double* ) malloc(N * sizeof(double));
        streamdata = readObsFile(data.obsFPath, N);

        printf("number of observations: %d\n", N);
    } else if (p == 2) {
        streamdata1 = readObsFile(data.obsFPath, N / p);
        streamdata2 = readObsFile(data.obsFPath2, N / p);

        streamdata = (double* ) malloc(N * sizeof(double));
        int j = 0;
        int i = 0;
        for (i = 0; i < N / p; i++) {
            streamdata[j++] = streamdata1[i];
            streamdata[j++] = streamdata2[i];
        }
        printf("number of observations: %d\n", j);
    } else if (p == 3) {
        streamdata1 = readObsFile(data.obsFPath, N / p);
        streamdata2 = readObsFile(data.obsFPath2, N / p);
        streamdata3 = readObsFile(data.obsFPath3, N / p);

        streamdata = (double* ) malloc(N * sizeof(double));
        int j = 0;
        int i = 0;
        for (i = 0; i < N / p; i++) {
            streamdata[j++] = streamdata1[i];
            streamdata[j++] = streamdata2[i];
            streamdata[j++] = streamdata3[i];
        }
        printf("number of observations: %d\n", j);
    }

    //	locations_obs_zsort_inplace(N, locations, streamdata);
    if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "diag_approx") == 0)
        EXAGEOSTAT_MLE_sdzcpy(&data, streamdata);

    if (log == 1 && test == 1)
        init_log(&data);

    START_TIMING(data.total_exec_time);
    nlopt_set_max_objective(opt, MLE_alg, (void *) &data);
    nlopt_optimize(opt, starting_theta, &opt_f);
    STOP_TIMING(data.total_exec_time);

    if (strcmp(data.actualZLocFPath, "") != 0) {
        printf("%s ========\n", data.actualZLocFPath);
        nZmiss = countlines(data.actualZLocFPath);
        missing_locations = readLocsFile(data.actualZLocFPath, N);
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
            printf("nZobs: %d, nZmiss:%d\n", nZobs, nZmiss);
            nZobs = N;
            Zobs = (double* ) malloc(p * nZobs * sizeof(double));
            Zactual = (double* ) malloc(p * nZmiss * sizeof(double));
            Zmiss = (double* ) malloc(p * nZmiss * sizeof(double));
            if (p == 2) {
                double* Zactual1 = (double* ) malloc(nZmiss * sizeof(double));
                double* Zactual2 = (double* ) malloc(nZmiss * sizeof(double));
                Zactual1 = readObsFile(data.actualZFPath, nZmiss);
                Zactual2 = readObsFile(data.actualZFPath2, nZmiss);
                int j = 0;
                int i = 0;
                for (i = 0; i < nZmiss; i++) {
                    Zactual[j++] = Zactual1[i];
                    Zactual[j++] = Zactual2[i];
                }
                printf("number of observations: %d\n", j);
                free(Zactual1);
                free(Zactual2);
            } else if (p == 3) {
                double* Zactual1 = (double* ) malloc(nZmiss * sizeof(double));
                double* Zactual2 = (double* ) malloc(nZmiss * sizeof(double));
                double* Zactual3 = (double* ) malloc(nZmiss * sizeof(double));
                Zactual1 = readObsFile(data.actualZFPath, nZmiss);
                Zactual2 = readObsFile(data.actualZFPath2, nZmiss);
                Zactual3 = readObsFile(data.actualZFPath3, nZmiss);

                int j = 0;
                int i = 0;
                for (i = 0; i < nZmiss; i++) {
                    Zactual[j++] = Zactual1[i];
                    Zactual[j++] = Zactual2[i];
                    Zactual[j++] = Zactual3[i];
                }
                printf("number of observations: %d\n", j);
                free(Zactual1);
                free(Zactual2);
                free(Zactual3);
            } else if (p == 1)
                Zactual = readObsFile(data.actualZFPath, nZmiss);
            MLE_get_zobs(&data, Zobs, N);
            data.lmiss = *missing_locations;
            data.lobs = *locations;

            printf("NZobs=%d, nZmiss=%d, N=%d,\n", nZobs, nZmiss, N);

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

        //***************
        if (data.mspe == 1) {
            //memory allocation
            Zobs = (double* ) malloc(p * nZobs * sizeof(double));
            Zactual = (double* ) malloc(p * nZmiss * sizeof(double));
            Zmiss = (double* ) malloc(p * nZmiss * sizeof(double));
            if (strcmp(data.computation, "exact") == 0 || strcmp(data.computation, "diag_approx") == 0)
                prediction_init(&data, nZmiss, nZobs, dts, p_grid, q_grid, 1);
#if defined(EXAGEOSTAT_USE_HICMA)
            else if (strcmp(data.computation, "lr_approx") == 0)
                prediction_init(&data, nZmiss, nZobs, lts, p_grid, q_grid, 1);
#endif

            int j = 0;
            for (j = 0; j < pred_samples; j++) {
                printf("nZobs = %d\n", p * nZobs);
                if (strcmp(data.actualZLocFPath, "") == 0) {
                    if (p == 2)
                        pick_random_points2(&data, Zobs, Zactual, nZmiss, nZobs, N);
                    else if (p == 3)
                        pick_random_points3(&data, Zobs, Zactual, nZmiss, nZobs, N);
                    else
                        pick_random_points(&data, Zobs, Zactual, nZmiss, nZobs, N);
                } else {
                    //	pred_samples=1;
                    pred_samples = 1;
                    if (p == 2) {
                        double* Zactual1 = (double* ) malloc(nZmiss * sizeof(double));
                        double* Zactual2 = (double* ) malloc(nZmiss * sizeof(double));
                        Zactual1 = readObsFile(data.actualZFPath, nZmiss);
                        Zactual2 = readObsFile(data.actualZFPath2, nZmiss);
                        int j = 0;
                        int i = 0;
                        for (i = 0; i < nZmiss; i++) {
                            Zactual[j++] = Zactual1[i];
                            Zactual[j++] = Zactual2[i];
                        }
                        printf("number of observations: %d\n", j);
                        free(Zactual1);
                        free(Zactual2);
                    } else if (p == 3) {
                        double* Zactual1 = (double* ) malloc(nZmiss * sizeof(double));
                        double* Zactual2 = (double* ) malloc(nZmiss * sizeof(double));
                        double* Zactual3 = (double* ) malloc(nZmiss * sizeof(double));
                        Zactual1 = readObsFile(data.actualZFPath, nZmiss);
                        Zactual2 = readObsFile(data.actualZFPath2, nZmiss);
                        Zactual3 = readObsFile(data.actualZFPath3, nZmiss);
                        int j = 0;
                        int i = 0;
                        for (i = 0; i < nZmiss; i++) {
                            Zactual[j++] = Zactual1[i];
                            Zactual[j++] = Zactual2[i];
                            Zactual[j++] = Zactual3[i];
                        }
                        printf("number of observations: %d\n", j);
                        free(Zactual1);
                        free(Zactual2);
                        free(Zactual3);
                    } else if (p == 1)
                        Zactual = readObsFile(data.actualZFPath, nZmiss);
                    MLE_get_zobs(&data, Zobs, N);
                    data.lmiss = *missing_locations;
                    data.lobs = *locations;
                }

                START_TIMING(pred_time);
                //generate_interior_points(&data, Zobs, NULL, nZmiss, nZobs, N);
                if (strcmp(data.computation, "exact") == 0)
                    prediction_error = EXAGEOSTAT_dmle_Predict_Tile(&data, starting_theta, p * nZmiss, p * nZobs, Zobs,
                                                                   Zactual, Zmiss, N);
                else if (strcmp(data.computation, "diag_approx") == 0)
                    prediction_error = EXAGEOSTAT_dmle_diag_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs,
                                                                        Zactual, Zmiss, N);
#if defined(EXAGEOSTAT_USE_HICMA)
                else if (strcmp(data.computation, "lr_approx") == 0) {
                    data.hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_POINT;
                    prediction_error = EXAGEOSTAT_TLR_dmle_Predict_Tile(&data, starting_theta, nZmiss, nZobs, Zobs, Zactual,
                                                               Zmiss, N, lts);
                }
#endif
                STOP_TIMING(pred_time);
                int index = 0;
                //for (index=0; index< nZmiss; index++)
                //	printf ("(%f, %f)\n ", Zactual[index], Zmiss[index]);

                fprintf(stderr, "Prediction Error: %e \n", prediction_error);
                avg_pred_value += prediction_error;
                avg_pred_value1 += data.mserror1;
                avg_pred_value2 += data.mserror2;
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

    print_result(&data, starting_theta, N, zvecs, ncores, lts, test, initial_theta, data.computation, p_grid, q_grid,
                 data.final_loglik, prediction_error);

    if (log == 1 && test == 1)
        finalize_log(&data);

    nlopt_destroy(opt);
    MLE_Finalize(&data);

    if (arguments.profile == 1) {
        starpu_fxt_stop_profiling();
        RUNTIME_profiling_display_efficiency();
    }

    return 0;
}
