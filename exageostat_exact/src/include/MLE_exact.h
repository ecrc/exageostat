/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_exact.h
 *
 * ExaGeoStat exact computation main functions header file.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#ifndef _MLE_EXACT_H_
#define _MLE_EXACT_H_

#include "MLE_misc.h"
#include "flat_file.h"

#if defined(EXAGEOSTAT_USE_NETCDF)

#include "nc_file.h"

#endif

#include "starpu_exageostat.h"


void EXAGEOSTAT_MLE_dzvg_Tile(MLE_data *data, double*Nrand, double*initial_theta,
                             int n, int dts, int log);

void EXAGEOSTAT_MLE_dzcpy(MLE_data *data, double*streamdata);

void EXAGEOSTAT_MLE_dzvg_Tile_Async(MLE_data *data, double*Nrand, double*initial_theta,
                                   int n, int dts, int log);

double EXAGEOSTAT_dmle_Tile(unsigned n, const double*theta, double*grad,
                           void *CHAM_data);

double EXAGEOSTAT_dmle_Tile_Async(unsigned n, const double*theta, double*grad,
                                 void *CHAM_data);

double EXAGEOSTAT_dmle_Predict_Tile(MLE_data *CHAM_data, double*theta, int nZmiss,
                                   int nZobs, double*Zobs, double*Zactual,
                                   double*Zmiss, int n);

void EXAGEOSTAT_dmle_Predict_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs,
                                     int dts, int p_grid, int q_grid,
                                     int mse_flag);

double EXAGEOSTAT_dmle_Predict_Tile_Async(MLE_data *CHAM_data, double*theta, int nZmiss,
                                         int nZobs, int n);

void EXAGEOSTAT_dmle_Call(MLE_data *data, int ncores, int gpus,
                         int dts, int p_grid, int q_grid,
                         int N, int nZobs, int nZmiss);

void EXAGEOSTAT_MLE_szcpy(MLE_data *data, double*streamdata);

void EXAGEOSTAT_psdpotrf(CHAM_enum uplo, CHAM_desc_t *A, int diag_thick,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void EXAGEOSTAT_dmle_mloe_mmom_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs,
                                       int dts, int p_grid, int q_grid);

void EXAGEOSTAT_dmle_mloe_mmom_Tile(MLE_data *CHAM_data, double*truth_theta, double*estimatedtheta,
                                   int nZmiss, int nZobs, int n);

void EXAGEOSTAT_dmle_mloe_mmom_Tile_Async(MLE_data *CHAM_data, double*truth_theta, double*estimatedtheta,
                                         int nZmiss, int nZobs, int n);

double*EXAGEOSTAT_Fisher_Tile(MLE_data *data, int N, double*initial_theta,
                              int dts, int p_grid, int q_grid);

#endif