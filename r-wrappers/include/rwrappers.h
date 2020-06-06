/**
 *
 * Copyright (c) 2017-2019  King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2019-07-30
 *
 **/
#ifndef _RWRAPPERS_H_
#define _RWRAPPERS_H_
#include "../../src/include/MLE.h"
/** ****************************************************************************
 *  Structure for reading real datasets through STARS-H
 **/
/*
   typedef struct {
   double x;                               ///< Locations x-axis.
   double y;                               ///< Locations y-axis.
   double z;                               ///< Measurements.
   } double;
   typedef struct {
   double theta;                               ///< Locations x-axis.
   double loglik;                               ///< Locations y-axis.
   int    iters;                               ///< Measurements.
   double time_in_secs;
   } double;
   */
void  gen_z_givenlocs_exact(double *x, int *xlen, double *y,
        int *ylen,  double *theta1,
        double *theta2, double *theta3,
        int *dmetric, int *n,
        int *ncores, int *gpus, int *ts,
        int *p_grid, int *q_grid, int *veclen,
        double *globalvec);

void  gen_z_exact(double *theta1, double *theta2,
        double *theta3, int *dmetric, int *n, int *seed,
        int *ncores,  int *gpus,
        int *ts, int *p_grid, int *q_grid,
        int *veclen,  double *globalvec);


void  mle_exact_non_stat(double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *dmetric, int *n, double *opt_tol,
        int *opt_max_iters, int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout);


void  mle_exact(double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *dmetric, int *n, double *opt_tol,
        int *opt_max_iters, int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout);

void  mle_tlr(  double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *tlr_acc, int *tlr_maxrank,
        int *dmetric, int *n,  double *opt_tol,
        int *opt_max_iters, int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout);

void  mle_dst(  double *x, int *xlen, double *y,
        int *ylen, double *z, int *zlen,
        double *clb, int *clblen, double *cub,
        int *cublen, int *dst_thick, int *dmetric,
        int *n, double *opt_tol, int *opt_max_iters,
        int *ncores, int *gpus,
        int *ts, int *p_grid, int *q_grid,
        double *globalthetaout);

void  rexageostat_init(int *ncores, int *gpus, int *ts);
void  rexageostat_finalize();

#endif
