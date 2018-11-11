/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
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
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _RWRAPPERS_H_
#define _RWRAPPERS_H_
#include "../../src/include/MLE.h"

void  gen_z_givenlocs_exact(int *n, int *ncores,  int *gpus,  int *ts,  int *p_grid, int *q_grid, double *x, int *xlen, double *y, int *ylen,  double *theta1, double *theta2, double *theta3, int *dmetric, int *veclen,  double *globalvec);

void  gen_z_exact(int *n, int *ncores,  int *gpus,  int *ts,  int *p_grid, int *q_grid,  double *theta1, double *theta2, double *theta3, int *dmetric, int *seed, int *veclen,  double *globalvec);

void  mle_exact(int *n,  int *ncores, int *gpus, int *ts, int *p_grid, int *q_grid,  double *x, int *xlen, double *y, int *ylen, double *z, int *zlen, double *clb, int *clblen, double *cub, int *cublen, int *dmetric, double *opt_tol, int *opt_max_iters, double *globalthetaout);

void  mle_tlr(int *n,  int *ncores, int *gpus, int *ts, int *p_grid, int *q_grid,  double *x, int *xlen, double *y, int *ylen, double *z, int *zlen, double *clb, int *clblen, double *cub, int *cublen, int *tlr_acc, int *tlr_maxrank, int *dmetric,  double *opt_tol, int *opt_max_iters, double *globalthetaout);

void  mle_dst(int *n,  int *ncores, int *gpus, int *ts, int *p_grid, int *q_grid,  double *x, int *xlen, double *y, int *ylen, double *z, int *zlen, double *clb, int *clblen, double *cub, int *cublen, int *dst_thick, int *dmetric, double *opt_tol, int *opt_max_iters, double *globalthetaout);

void  rexageostat_init(int *ncores, int *gpus, int *ts);
void  rexageostat_finalize();

#endif
