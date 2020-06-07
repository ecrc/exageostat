/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#ifndef _MLE_H_
#define _MLE_H_
#include "MLE_exact.h"
#include "MLE_exact_s.h"
#include "MLE_approx.h"
#include "MLE_sdexact.h"
#if defined( EXAGEOSTAT_USE_HICMA )
#include "MLE_lr.h"
#include "hicma_constants.h"
#endif
/** ****************************************************************************
 *  Macro to print a warning message to the user about the propoer way to input his arguments
 **/
#define USAGE(args, details)  fprintf(stderr," Proper Usage is : ./examples/zgen_mle_test "args" \n"   details);

void MLE_zvg(MLE_data *data,  double *Nrand, double *initial_theta,
        int n, int ts, int log,
        int p_grid, int q_grid);
//void MLE_zvr(MLE_data *data, int n, char *format);
//
double MLE_alg(unsigned n, const double *theta, double *grad,
        void *data);

void exageostat_init(int *ncores, int *gpus, int *dts,
        int *lts);

void exageostat_finalize();

void prediction_init(MLE_data *data, int nZmiss, int nZobs,
        int ts, int p_grid, int q_grid,
        int mse_flag);

void prediction_finalize(MLE_data *data);

void MLE_Finalize(MLE_data *data);

void MLE_get_zobs(MLE_data *data, double *z, int n);

void MLE_sdregister(MLE_data *data);

void mloe_mmom_init(MLE_data *data, int nZmiss, int nZobs,
        int ts, int p_grid, int q_grid);

void MLOE_MMOM_Finalize(MLE_data *data);
#endif
