/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file synthetic_sdmle_rwrapper_test
 *
 * ExaGeoStat main functions.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "examples.h"
#include "../src/include/MLE.h"
#include "../r-wrappers/include/rwrappers.h"

int main(int argc, char** argv) {

    //initialization
    char* theta;  //for testing case
    int n, dts, lts;//, log, verbose;
    int i = 0, seed = 0;
    char* dm;
    char* computation, *clb, *cub;
    char* kernel_fun;
    int dm_int = 0, hicma_acc = 0, hicma_maxrank = 0, diag_thick = 0, num_params = 0;
    int globalveclen;
    int gpus = 0, p_grid, q_grid, ncores;
    arguments arguments;
    double* vecs_out = NULL, *theta_out = NULL;
    double* initial_theta = NULL, *lb = NULL, *ub = NULL;
    double opt_tol;
    int opt_max_iters = 0;
    //Arguments default values
    set_args_default(&arguments);
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    //read inputs
    n = atoi(arguments.N);
    ncores = atoi(arguments.ncores);
    gpus = atoi(arguments.gpus);
    p_grid = atoi(arguments.p);
    q_grid = atoi(arguments.q);
    diag_thick = atoi(arguments.diag_thick);
    dts = atoi(arguments.dts);
    lts = atoi(arguments.lts);
    hicma_acc = atoi(arguments.acc);
    hicma_maxrank = atoi(arguments.maxrank);
    dm = arguments.dm;
    computation = arguments.computation; //approx or exact
    theta = arguments.ikernel;
    clb = arguments.olb;
    cub = arguments.oub;
    opt_tol = pow(10, -1.0 * atoi(arguments.opt_tol));
    opt_max_iters = atoi(arguments.opt_max_iters);
    kernel_fun = arguments.kernel_fun;
    dm_int = strcmp(dm, "ed") == 0 ? 0 : 1;

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
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        num_params = 6;
    else if (strcmp(kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
        num_params = 6;
    else if (strcmp(kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        num_params = 7;
    else if (strcmp(kernel_fun, "univariate_matern_non_stat") == 0)
        num_params = 8;
    else {
        fprintf(stderr, "Choosen kernel is not exist!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    globalveclen = num_params * n;
    //Memory allocation
    theta_out = (double* ) malloc((num_params + 3) * sizeof(double));
    vecs_out = (double* ) malloc(globalveclen * sizeof(double));
    initial_theta = (double* ) malloc(num_params * sizeof(double));
    lb = (double* ) malloc(num_params * sizeof(double));
    ub = (double* ) malloc(num_params * sizeof(double));
    theta_parser2(initial_theta, theta, num_params);
    theta_parser2(lb, clb, num_params);
    theta_parser2(ub, cub, num_params);

    rexageostat_init(&ncores, &gpus, &dts, &lts);

    int kernel_len = strlen("univariate_matern_stationary");

    gen_z_exact("univariate_matern_stationary", initial_theta,
                &num_params,
                &dm_int, &n, &seed, &ncores,
                &gpus, &dts, &p_grid,
                &q_grid, &globalveclen, vecs_out);

    if (strcmp(computation, "exact") == 0)
        mle_mp(vecs_out, NULL, &vecs_out[n],
               NULL, &vecs_out[2 * n], NULL,
               lb, &num_params, ub,
               &num_params, &diag_thick, "univariate_matern_stationary",
               &dm_int, &n,
               &opt_tol, &opt_max_iters, &ncores,
               &gpus, &lts, &p_grid, &q_grid,
               theta_out);


    rexageostat_finalize();
    for (int i = 0; i < num_params; i++)
        printf("%f -", theta_out[i]);
    printf("\n");

    printf("time_per_iter: %f secs\n", theta_out[num_params]);
    printf("total_time: %f secs\n", theta_out[num_params + 1]);
    printf("no_iters: %f \n", theta_out[num_params + 2]);
    //free memory
    free(theta_out);
    free(vecs_out);
    free(ub);
    free(lb);
    return 0;
}