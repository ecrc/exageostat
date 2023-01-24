/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE.h
 *
 * Header file of ExaGeoStat main functions.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2021-07-28
 *
 **/

#ifndef _EXAGEOSTAT_EXAMPLES_H_
#define _EXAGEOSTAT_EXAMPLES_H_

#include <stdlib.h>
#include <argp.h>
#include "MLE_misc.h"

/** ****************************************************************************
 *  EXAGEOSTAT arguments uniquely identifies a set of input arguments
 **/
typedef struct {
    int test;                   ///< The -t flag -- use testing mode (Synthetic Dataset Generator).
    int check;                  ///< The -c flag -- use check mode (in the case of approximation only).
    char *zvecs;                ///< The number of Z vectors to be tested.
    int verbose;                ///< The -v flag -- verbose mode.
    char *ncores;               ///< Number of CPU cores.
    char *gpus;                 ///< Number of GPUs.
    char *N;                    ///< Problem size -- Only in the case of testing mode.
    char *p;                    ///< P in distributed grid.
    char *q;                    ///< q in distributed grid.
    char *lts;                  ///< HiCMA tile size.
    char *dts;                  ///< Chameleon tile size.
    char *kernel;               ///< Target theta vector output (ex., ?:?:?).
    char *ikernel;              ///< Initial theta vector -- Only in the case of testing mode (ex., 1:0.1:0.5).
    char *olb;                  ///< Optimizer lower bounds vector (ex., 0.1:0.01:0.01).
    char *oub;                  ///< Optimizer upper bounds vector (ex., 5:5:5).
    char *computation;          ///< Approx or exact.
    int async;                  ///< 0--> tile  1--> tile_async.The --async flag.
    char *locs_file;            ///< Locations file path -- in the case of real dataset (real mode).
    char *time_file;            ///< Time file path -- in the case of real dataset (real mode -- space-time kernel).
    char *obs_dir;              ///< Observations file path --  in the case of real dataset (real mode).
    char *obs_dir2;             ///< Observations file path2 (bivariate case) --  in the case of real dataset (real mode).
    char *obs_dir3;             ///< Observations file path3 (bivariate case) --  in the case of real dataset (real mode).
    char *actualZ_file;         ///< Actual observations file path -- in the case of prediction.
    char *actualZ_file2;        ///< Actual observations file path -- in the case of prediction.
    char *actualZ_file3;        ///< Actual observations file path -- in the case of prediction.
    char *actualZloc_file;      ///< Actual locations file path -- in the case of prediction.
    char *actualtime_file;      ///< Actual time file path -- in the case of prediction (space-time kernel).
    char *predict;              ///< Number of missing values  -- in the case of testing mode.
    char *dm;                   ///< Distance metric to be used ed->Euclidian Distance -- gcd->Greate Circle Distance.
    char *diag_thick;           ///< The thick of used diagonal in the case of diagonal approximation approch.
    int log;                    ///< Generate log files -- 0->do not store generated data, 1-->store generated data.
    char *maxrank;              ///< Max Rank in the case of LR-HiCMA approx.
    char *acc;                  ///< Accuracy in the case of LR-HiCMA approx.
    int profile;                ///< profiling the performance of exageostat using FxT.
    char *opt_tol;              ///< The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
    char *opt_max_iters;        ///< Maximum number of mle iterations.
    int ooc;                    ///< Support Out-Of-Core (OOC) -- 0->do not support, 1-->support.
    char *kernel_fun;           ///< Stationary_matern, or non_stationary_matern.
    int mloe_mmom;              ///< Use MLOE and MMOM
    int mloe_mmom_async;        ///< Use MLOE and MMOM Async
    int mspe;                   ///>compute mspe.
    int idw;                    ///<IDW prediction.
    char *checkpoint_file;      ///< checkpoint file path.
    char *recovery_file;        ///< Recovery file path.
    char *dim;                  ///< 2D or 3D
    char *time_slots;           ///< spatio-temoral time slots.
    int fisher;                 ///< fisher matrix
} arguments;

void check_args(arguments *arg_values);

void set_args_default(arguments *arg_values);

void init(int *test, int *N, int *ncores,
          int *gpus, int *p_grid, int *q_grid,
          int *zvecs, int *dts, int *lts,
          int *nZmiss, int *log, double *initial_theta,
          double *starting_theta, double *target_theta, double *lb,
          double *ub, MLE_data *data, arguments *arguments);

static struct argp_option options[] =
        {
                {"test",            'a', 0,                         0, "Execute in test mode"},
                {"check",           'b', 0,                         0, "Produce check output"},
                {"zvecs",           'c', "ZVECS",                   0, "number of Z vectors to be tested"},
                {"verbose",         'd', 0,                         0, "Produce verbose output"},
                {"ncores",          'e', "NCORES",                  0, "Number of cores"},
                {"gpus",            'f', "GPUS",                    0, "Number of gpus"},
                {"p",               'g', "P",                       0, "p in distributed system"},
                {"q",               'h', "Q",                       0, "q in distributed system"},
                {"N",               'i', "MATRIX_SIZE",             0, "Synthetic Matrix Size"},
                {"lts",             'j', "HICMA_TILE_SIZE",         0, "Number of tiles in TLR"},
                {"dts",             'k', "Chameleon_TILE_SIZE",     0, "Number of Tiles in dense"},
                {"kernel",          'l', "KERNEL",                  0, "Computation model"},
                {"ikernel",         'm', "IKERNEL",                 0, "Initial theta(s) used for testing case generation"},
                {"olb",             'n', "LB",                      0, "Optimizer Lower Bounds"},
                {"oub",             'o', "UB",                      0, "Optimizer Upper Bounds"},
                {"computation",     'p', "COMPUTATION",             0, "Exact or Approx"},
                {"async",           'q', 0,                         0, "Asynchronous"},
                {"locs_file",       'r', "LOCATIONS_FILE",          0, "Read Locations from this Location File"},
                {"obs_dir",         's', "OBSERVATIONS_DIRECTORY",  0, "Read Observations from this directory path"},
                {"obs_dir2",        't', "OBSERVATIONS_DIRECTORY2", 0, "Read Observations from this directory path"},
                {"actualZ_file",    'u', "ACTUALZ_FILE",            0, "Read actual Z from this observation file"},
                {"actualZ_file2",   'v', "ACTUALZ2_FILE",           0, "Read actual Z from this observation file2"},
                {"actualZloc_file", 'w', "ACTUALZLOC_FILE",         0, "Read actual Z locations from this location file"},
                {"predict",         'x', "PREDICT",                 0, "Number of Missing Values"},
                {"dm",              'y', "DISTANCE_METRIC",         0, "Distance Metric"},
                {"diag_thick",      'z', "DIAG_THICK",              0, "Diagonal Thick"},
                {"log",             'A', 0,                         0, "Store Generated Data (Test mode only)"},
                {"maxrank",         'B', "MAXRANK",                 0, "HiCMA Max RANK"},
                {"acc",             'C', "ACC",                     0, "HiCMA Accuracy"},
                {"profile",         'D', 0,                         0, "Performance profiling"},
                {"opt_tol",         'E', "OPTIMIZATION_TOLERANCE",  0, "Optimization tolerance"},
                {"opt_iters",       'F', "OPTIMIZATION_MAX_ITERS",  0, "Optimization maximum iterations"},
                {"ooc",             'G', 0,                         0, "Support Out-Of-Core (OOC) execution"},
                {"kernel_fun",      'H', "Core Kernel",             0, "stationary_matern or nonstationary_matern or bivariate_kernel"},
                {"mloe_mmom",       'I', 0,                         0, "use_mloe_mmom"},
                {"mloe_mmom_async", 'J', 0,                         0, "use_mloe_mmom_async"},
                {"mspe",            'K', 0,                         0, "use_mspe"},
                {"checkpoint_file", 'L', "CHECKPOINT_FILE",         0, "Checkpoint parameters to"},
                {"recovery_file",   'M', "RECOVERY_FILE",           0, "Recover parameters from"},
                {"dim",             'N', "DIMENSION",               0, "Dimension"},
                {"time_slots",      'O', "TIME_SLOT",               0, "time slots"},
                {"idw",             'P', 0,                         0, "inverse distance weighted (IDW)"},
                {"time_file",       'Q', "TIME_FILE",               0, "Read time slots from this time File"},
                {"actualtime_file", 'R', "ACTUALZTIME_FILE",        0, "Read actual time slots from this time file"},
                {"obs_dir3",        'S', "OBSERVATIONS_DIRECTORY3", 0, "Read Observations from this directory path"},
                {"actualZ_file3",   'T', "ACTUALZ3_FILE",           0, "Read actual Z from this observation file3"},
                {"fisher",          'U', 0,                         0, "Compute Fisher matrix"},
                {0}
        };

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    arguments *arguments = state->input;

    switch (key) {
        case 'a':
            arguments->test = 1;
            break;
        case 'b':
            arguments->check = 1;
            break;
        case 'c':
            arguments->zvecs = arg;
            break;
        case 'd':
            arguments->verbose = 1;
            break;
        case 'e':
            arguments->ncores = arg;  //non-optional;
            break;
        case 'f':
            arguments->gpus = arg;  //non-optional;
            break;
        case 'g':
            arguments->p = arg;  //non-optional;
            break;
        case 'h':
            arguments->q = arg;
            break;
        case 'i':
            arguments->N = arg;
            break;
        case 'j':
            arguments->lts = arg;  //non-optional
            break;
        case 'k':
            arguments->dts = arg;  //non-optional
            break;
        case 'l':
            arguments->kernel = arg;
            break;
        case 'm':
            arguments->ikernel = arg;
            break;
        case 'n':
            arguments->olb = arg;
            break;
        case 'o':
            arguments->oub = arg;
            break;
        case 'p':
            arguments->computation = arg;
            break;
        case 'q':
            arguments->async = 1;
            break;
        case 'r':
            arguments->locs_file = arg;
            break;
        case 's':
            arguments->obs_dir = arg;
            break;
        case 't':
            arguments->obs_dir2 = arg;
            break;
        case 'u':
            arguments->actualZ_file = arg;
            break;
        case 'v':
            arguments->actualZ_file2 = arg;
            break;
        case 'w':
            arguments->actualZloc_file = arg;
            break;
        case 'x':
            arguments->predict = arg;
            break;
        case 'y':
            arguments->dm = arg;
            break;
        case 'z':
            arguments->diag_thick = arg;
            break;
        case 'A':
            arguments->log = 1;
            break;
        case 'B':
            arguments->maxrank = arg;
            break;
        case 'C':
            arguments->acc = arg;
            break;
        case 'D':
            arguments->profile = 1;
            break;
        case 'E':
            arguments->opt_tol = arg;
            break;
        case 'F':
            arguments->opt_max_iters = arg;
            break;
        case 'G':
            arguments->ooc = 1;
            break;
        case 'H':
            arguments->kernel_fun = arg;
            break;
        case 'I':
            arguments->mloe_mmom = 1;
            break;
        case 'J':
            arguments->mloe_mmom_async = 1;
            break;
        case 'K':
            arguments->mspe = 1;
            break;
        case 'L':
            arguments->checkpoint_file = arg;
            break;
        case 'M':
            arguments->recovery_file = arg;
            break;
        case 'N':
            arguments->dim = arg;
            break;
        case 'O':
            arguments->time_slots = arg;
            break;
        case 'P':
            arguments->idw = 1;
            break;
        case 'Q':
            arguments->time_file = arg;
            break;
        case 'R':
            arguments->actualtime_file = arg;
            break;
        case 'S':
            arguments->obs_dir3 = arg;
            break;
        case 'T':
            arguments->actualZ_file3 = arg;
            break;
        case 'U':
            arguments->fisher = 1;
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static char args_doc[] = "";

static char doc[] =
        "ExaGeoStat -- A unified geospatial statistic framework to evaluate Maximum Likelihood function using both real and synthetic dataset on CHAMELEON (Dense/DST/mixed-precision) - HiCMA (TLR approximation) - )";

static struct argp argp = {options, parse_opt, args_doc, doc};

#endif