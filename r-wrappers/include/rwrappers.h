/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file rwrappers.h
 *
 * ExaGeoStat R-wrapper functions header.
 *
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-14
 *
 **/
#ifndef _RWRAPPERS_H_
#define _RWRAPPERS_H_
#include "../../src/include/MLE.h"

void  rexageostat_gen_z(int *n, int *ncores,  int *gpus,  int *ts,  int *p_grid, int *q_grid,  double *theta1, double *theta2, double *theta3,  int *computation, int *dmetric, int *veclen,  double *globalvec);
void  rexageostat_likelihood(int *n,  int *ncores, int *gpus, int *ts, int *p_grid, int *q_grid,  double *x, int *xlen, double *y, int *ylen, double *z, int *zlen, double *clb, int *clblen, double *cub, int *cublen,  int *computation, int *dmetric, double *opt_tol, int *opt_max_iters, double *globalthetaout);
void  rexageostat_init(int *ncores, int *gpus, int *ts);
void  rexageostat_finalize();

#endif
