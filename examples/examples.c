/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file examples.c
 *
 *
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "examples.h"

output results;

void set_args_default(arguments *arg_values)
//! set default values for input
/*!  arguments
 * @param[in] arg_values: user arguments
 * */
{
    arg_values->test = 0;
    arg_values->check = 0;
    arg_values->verbose = 0;
    arg_values->zvecs = "1";
    arg_values->computation = "exact";
    arg_values->async = 0;
    arg_values->kernel = "";
    arg_values->ikernel = "";
    arg_values->ncores = "1";
    arg_values->gpus = "0";
    arg_values->p = "1";
    arg_values->q = "1";
    arg_values->N = "0";
    arg_values->lts = "0";
    arg_values->dts = "0";
    arg_values->locs_file = "";
    arg_values->obs_dir = "";
    arg_values->obs_dir2 = "";
    arg_values->actualZ_file = "";
    arg_values->actualZloc_file = "";
    arg_values->predict = "0";
    arg_values->dm = "ed";
    arg_values->diag_thick = "1";
    arg_values->log = 0;
    arg_values->maxrank = "0";
    arg_values->acc = "0";
    arg_values->profile = 0;
    arg_values->opt_tol = "5";
    arg_values->opt_max_iters = "-1";
    arg_values->ooc = 0;
    arg_values->kernel_fun = "univariate_matern_stationary";
    arg_values->mloe_mmom = 0;
    arg_values->mloe_mmom_async = 0;
    arg_values->mspe = 0;
    arg_values->test = 0;
    arg_values->check = 0;
    arg_values->verbose = 0;
    arg_values->zvecs = "1";
    arg_values->computation = "exact";
    arg_values->computation = "matern";
    arg_values->async = 0;
    arg_values->kernel = "";
    arg_values->ikernel = "";
    arg_values->ncores = "1";
    arg_values->gpus = "0";
    arg_values->p = "1";
    arg_values->q = "1";
    arg_values->N = "0";
    arg_values->lts = "0";
    arg_values->dts = "0";
    arg_values->locs_file = "";
    arg_values->obs_dir = "";
    arg_values->actualZ_file = "";
    arg_values->actualZ_file2 = "";
    arg_values->actualZloc_file = "";
    arg_values->predict = "0";
    arg_values->dm = "ed";
    arg_values->diag_thick = "1";
    arg_values->log = 0;
    arg_values->maxrank = "0";
    arg_values->acc = "0";
    arg_values->profile = 0;
    arg_values->opt_tol = "5";
    arg_values->opt_max_iters = "-1";
    arg_values->ooc = 0;
    arg_values->recovery_file = "";
    arg_values->checkpoint_file = "";
    arg_values->dim = "2d";
    arg_values->time_slots = "0";
    arg_values->idw = 0;
    arg_values->seed = "0";
}

void check_args(arguments *arg_values) {
    //! check  values for input
    /*!  arguments
     * @param[in] arg_values: user arguments
     * */
    if (arg_values->test == 0) {
        if (strcmp(arg_values->locs_file, "") == 0 || strcmp(arg_values->obs_dir, "") == 0) {
            fprintf(stdout,
                    "(ExaGeoStat Error MSG): Real running mode requires both locs_file and obs_dir path.... \nExaGeoStat terminated \n\n");
            exit(EXIT_FAILURE);
        }
    } else {
        if (strcmp(arg_values->locs_file, "") != 0 || strcmp(arg_values->obs_dir, "") != 0)
            fprintf(stdout,
                    "(ExaGeoStat Warning MSG): Test running mode does not require locs_file and obs_dir paths, any will be ignored-> continue running...\n");
    }
}


void
init(int *test, int *N, int *ncores, int *gpus, int *p_grid, int *q_grid, int *zvecs, int *dts, int *lts, int *nZmiss,
     int *log, double* initial_theta, double* starting_theta, double* target_theta, double* lb, double* ub,
     MLE_data *data, arguments *arguments)
//! initialize exageostat by setting several
/*! variables from argument
 * Returns MLE_data struct.
 * @param[in] arguments: command line arguments.
 * @param[out] test: if test=1 ->test running mode  if test=0->real running mode.
 * @param[out] N: number of spatial locations (problem size).
 * @param[out] ncores: number of CPU computing units.
 * @param[out] gpus: number of GPU computing units.
 * @param[out] p_grid: p_grid in the case of distributed system.
 * @param[out] q_grid: q_grid in the case of distributed system.
 * @param[out] zvecs: the number of Z vectors that should be generated in the case of test running mode.
 * @param[out] dts: dense tile size.
 * @param[out] lts: TLR tile size.
 * @param[out] nZmiss: number of unknown observation to be predicted in the case of test running mode.
 * @param[out] log: determine if log files should be generated or not (log files stored on disk)
 * @param[out] initial_theta: initial_theta Vector with three parameter (Variance, Range, Smoothness)
 that is used to to generate the Covariance Matrix and initial Z vector.
 * @param[out] starting_theta: theta Vector with three parameter (Variance, Range, Smoothness)
 that is used to to generate the Covariance Matrix of the first MLE iteration.
 * @param[out] target_theta: target theta Vector with three parameter (Variance, Range, Smoothness) unknown theta parameter should be shown as '?'
 * @param[out] lb: optimization lower bounds vector ( lb_1, lb_2, lb_3).
 * @param[out] ub: optimization upper bounds vector ( ub_1, ub_2, ub_3).
 * @param[out] data: MLE_data struct with different MLE inputs.
 * @param[out] arguments: command line arguments.
 * */
{
    init_data_values(data);
    int i = 0;
    int num_params = 0;
    *test = arguments->test;
    *ncores = atoi(arguments->ncores);
    *gpus = atoi(arguments->gpus);
    *p_grid = atoi(arguments->p);
    *q_grid = atoi(arguments->q);
    *N = atoi(arguments->N);
    *zvecs = atoi(arguments->zvecs);
    *dts = atoi(arguments->dts);
    *lts = atoi(arguments->lts);
    *nZmiss = atoi(arguments->predict);
    *log = arguments->log;
    data->computation = arguments->computation;// exact or approx.
    data->async = arguments->async; // 0-->tile  1-->tile_async.
    data->locsFPath = arguments->locs_file;
    data->obsFPath = arguments->obs_dir;
    data->obsFPath2 = arguments->obs_dir2;
    data->obsFPath3 = arguments->obs_dir3;
    data->actualZFPath = arguments->actualZ_file;
    data->actualZFPath2 = arguments->actualZ_file2;
    data->actualZFPath3 = arguments->actualZ_file3;
    data->actualZLocFPath = arguments->actualZloc_file;
    data->timeFPath = arguments->time_file;
    data->actualTimeFPath = arguments->actualtime_file;
    data->dm = arguments->dm;
    data->diag_thick = atoi(arguments->diag_thick);
    data->log = arguments->log;
    data->hicma_maxrank = atoi(arguments->maxrank);
    data->hicma_acc = atof(arguments->acc);
    data->check = arguments->check;
    data->verbose = arguments->verbose;
    data->opt_tol = atoi(arguments->opt_tol);
    data->opt_max_iters = atoi(arguments->opt_max_iters);
    data->ooc = arguments->ooc;

    data->mloe_mmom = arguments->mloe_mmom;
    data->mloe_mmom_async = arguments->mloe_mmom_async;
    data->mspe = arguments->mspe;
    data->idw = arguments->idw;
    data->fisher = arguments->fisher;
    data->mloe_mmom = arguments->mloe_mmom;
    data->mloe_mmom_async = arguments->mloe_mmom_async;
    data->mspe = arguments->mspe;
    data->idw = arguments->idw;
    data->recovery_file = arguments->recovery_file;
    data->checkpoint_file = arguments->checkpoint_file;

    data->time_slots = atoi(arguments->time_slots);
    data->kernel_fun = arguments->kernel_fun;

    data->seed = atoi(arguments->seed);
    //should be improved
    location l;
    l.x = (double* ) malloc(sizeof(double));
    l.y = (double* ) malloc(sizeof(double));
    l.x[0] = 0.5;
    l.y[0] = 0.5;
    data->lm = l;
    //****************************************

    if (strcmp(data->kernel_fun, "univariate_matern_stationary") == 0 ||
        strcmp(data->kernel_fun, "univariate_pow_exp_stationary") == 0)
        num_params = 3;
    else if (strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        num_params = 4;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0)
        num_params = 9;
    else if (strcmp(data->kernel_fun, "bivariate_matern_flexible") == 0)
        num_params = 11;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        num_params = 7;
    else if (strcmp(data->kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        num_params = 10;
    else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0 ||
             strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0)
        num_params = 10;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stat") == 0)
        num_params = 8;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_gaussian") == 0 ||
             strcmp(data->kernel_fun, "univariate_exp_non_gaussian") == 0)
        num_params = 6;
    else {
        fprintf(stderr, "Choosen kernel is not exist(6)!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }

    theta_parser2(lb, arguments->olb, num_params);
    theta_parser2(ub, arguments->oub, num_params);


    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        lb[0] = ub[0] = 1;
        lb[1] = ub[1] = 1;
    } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {
        lb[0] = ub[0] = 1;
        lb[1] = ub[1] = 1;
        lb[2] = ub[2] = 1;
    }

    for (i = 0; i < num_params; i++)
        starting_theta[i] = lb[i];//(double)abs(lb[i])+(double)abs(ub[i])/2.0;

    theta_parser(initial_theta, target_theta, starting_theta, arguments->ikernel, arguments->kernel, lb, ub, *test,
                 num_params);

    //fill results value
    results.problem_size = *N;

    results.computation = malloc(strlen(data->computation) + 1);
    strcpy(results.computation, data->computation);

    results.kernel = malloc(strlen(data->kernel_fun) + 1);
    strcpy(results.kernel, data->kernel_fun);

    if (*test == 1)
        results.ds_type = "test";
    else
        results.ds_type = "real";

    if (data->precision == 0)
        results.precision = "double";
    else if (data->precision == 1)
        results.precision = "single";
    else if (data->precision == 2)
        results.precision = "double/single";

    results.dense_ts = *dts;

    results.lr_ts = *lts;

    results.ncores = *ncores;

    results.ngpus = *gpus;

    results.p = *p_grid;

    results.q = *q_grid;

    results.num_params = num_params;


    results.initial_theta = (double* ) malloc(num_params * sizeof(double));
    results.starting_theta = (double* ) malloc(num_params * sizeof(double));
    results.estimated_theta = (double* ) malloc(num_params * sizeof(double));

    for (i = 0; i < num_params; i++)
        results.initial_theta[i] = initial_theta[i];//(double)abs(lb[i])+(double)abs(ub[i])/2.0;

    for (i = 0; i < num_params; i++)
        results.starting_theta[i] = starting_theta[i];//(double)abs(lb[i])+(double)abs(ub[i])/2.0;

    results.lr_acc = data->hicma_maxrank;
    results.lr_maxrank = data->hicma_acc;


    results.lr_acc = data->hicma_maxrank;
    results.lr_maxrank = data->hicma_acc;

}
