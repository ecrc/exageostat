/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#ifndef _MLE_EXACT_H_
#define _MLE_EXACT_H_
#include "MLE_misc.h"
#include "starpu_exageostat.h"



void MORSE_MLE_dzvg_Tile (MLE_data *data,  double *Nrand, double *initial_theta, int n, int ts, int test, int log);
void MORSE_MLE_dzvg_Tile_Async (MLE_data *data,  double *Nrand, double *initial_theta, int n, int ts, int test, int log);
double MORSE_MLE_Tile(unsigned n, const double *theta, double *grad, void *MORSE_data);
double MORSE_MLE_Tile_Async(unsigned n, const double *theta, double *grad, void *MORSE_data);
double MORSE_MLE_Predict_Tile(MLE_data *MORSE_data, double *theta, int nZmiss, int nZobs, int n);
double MORSE_MLE_Predict_Tile_Async(MLE_data *MORSE_data, double *theta, int nZmiss, int nZobs, int n);
void MORSE_Call(MLE_data *data, int ncores,int gpus, int ts, int p_grid, int q_grid, int N, int nZobs, int nZmiss);


#endif
