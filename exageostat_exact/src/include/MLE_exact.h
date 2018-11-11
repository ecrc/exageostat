/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
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
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
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



void MORSE_MLE_dzvg_Tile (MLE_data *data,  double *Nrand, double *initial_theta, int n, int dts, int log);
void MORSE_MLE_zcpy( MLE_data *data, double *streamdata);
void MORSE_MLE_dzvg_Tile_Async (MLE_data *data,  double *Nrand, double *initial_theta, int n, int dts, int log);
double MORSE_MLE_Tile(unsigned n, const double *theta, double *grad, void *MORSE_data);
double MORSE_MLE_Tile_Async(unsigned n, const double *theta, double *grad, void *MORSE_data);
double MORSE_MLE_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n);
void MORSE_MLE_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid, int mse_flag);
double MORSE_MLE_Predict_Tile_Async(MLE_data *MORSE_data, double *theta, int nZmiss, int nZobs, int n);
void MORSE_Call(MLE_data *data, int ncores,int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss);


#endif
