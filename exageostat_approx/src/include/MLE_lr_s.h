/**
 *
 *  Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 *  All rights reserved.
 *  **
 * @file MLE.h
 *
 *
 *  ExaGeoStat is a software package provided by KAUST,
 *
 * @version 1.2.0
 * @author Sameh Abdulah
 * @date 2022-11-09
 * @generated d Fri Nov 22 15:11:13 2016
 *
 **/
#ifndef _MLE_LR_S_H_
#define _MLE_LR_S_H_

#include "MLE_misc.h"
#include "flat_file.h"

#if defined(EXAGEOSTAT_USE_NETCDF)

#include "nc_file.h"

#endif

#include "starpu_exageostat_approx.h"
#include "starsh-spatial.h"
#include "hicma_struct.h"
#include "include/hicma_z.h"
#include "misc/auxcompute_z.h"
#include "misc/auxdescutil.h"


//Generate Observations Vector (Z) for testing Maximum Likelihood function -- CHAMELEON-sync
int HICMA_MLE_szvg_Tile(MLE_data *data, double* Nrand, double* initial_theta,
                        int n, int dts, int log, int p_grid, int q_grid);

//Generate Observations Vector (Z) for testing Maximum Likelihood function -- CHAMELEON-async
int HICMA_MLE_szvg_Tile_Async(MLE_data *data, double* Nrand, double* initial_theta,
                              int n, int dts, int log, int p_grid, int q_grid);

//compute MLE function --syn
double HICMA_smle_Tile(unsigned n, const double* theta, double* grad, void *data);

//compute MLE function --async
double HICMA_smle_Tile_Async(unsigned n, const double* theta, double* grad, void *data);

void EXAGEOSTAT_TLR_dmle_Predict_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs,
                                 int lts, int p_grid, int q_grid, int mse_flag);

//Predict missing values base on a set of given values and covariance matrix
double HICMA_smle_Predict_Tile(MLE_data *CHAM_data, double* theta, int nZmiss,
                               int nZobs, double* Zobs, double* Zactual,
                               double* Zmiss, int n, int lts);

//init Hicma decriptors
void HICMA_smle_Call(MLE_data *data, int ncores, int ts,
                     int N, int nZobs, int nZmiss,
                     int gpus, int p_grid, int q_grid);

void EXAGEOSTAT_TLR_MLE_zcpy(MLE_data *data, double* streamdata);

#endif
