/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#ifndef _MLE_EXACT_S_H_
#define _MLE_EXACT_S_H_
#include "MLE_misc.h"
#include "flat_file.h"
#if defined(EXAGEOSTAT_USE_NETCDF)
#include "nc_file.h"
#endif
#include "starpu_exageostat.h"


void MORSE_MLE_szvg_Tile (MLE_data *data,  float *Nrand, double *initial_theta,
        int n, int dts, int log);
//void MORSE_MLE_zcpy( MLE_data *data, double *streamdata);

void MORSE_MLE_szvg_Tile_Async (MLE_data *data,  float *Nrand, double *initial_theta,
        int n, int dts, int log);

double MORSE_smle_Tile(unsigned n, const double *theta, double *grad,
        void *MORSE_data);

double MORSE_smle_Tile_Async(unsigned n, const double *theta, 
        double *grad, void *MORSE_data);

double MORSE_smle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss,
        int nZobs, double *Zobs, double *Zactual,
        double *Zmiss, int n);

void MORSE_smle_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs,
        int dts, int p_grid, int q_grid,
        int mse_flag);

double MORSE_smle_Predict_Tile_Async(MLE_data *MORSE_data, double * theta, 
        int nZmiss, int nZobs, double *Zobs,
        double *Zactual, double *Zmiss, int n);

void MORSE_smle_Call(MLE_data *data, int ncores,int gpus,
        int dts, int p_grid, int q_grid, 
        int N, int nZobs, int nZmiss);


#endif
