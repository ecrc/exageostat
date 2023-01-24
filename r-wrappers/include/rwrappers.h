/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
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
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#ifndef _RWRAPPERS_H_
#define _RWRAPPERS_H_

#include "../../src/include/MLE.h"

void gen_z_givenlocs_exact(double*x, int *xlen, double*y,
                           int *ylen, int *kernel, double* theta,
                           int *thetalen, int *dmetric, int *n,
                           int *ncores, int *gpus, int *ts,
                           int *p_grid, int *q_grid, int *veclen,
                           double* globalvec);

void gen_z_exact(int *kernel, double* theta, int *thetalen,
                 int *dmetric, int *n, int *seed,
                 int *ncores, int *gpus, int *ts,
                 int *p_grid, int *q_grid, int *veclen,
                 double* globalvec);

void mle_exact(double* x, int *xlen, double* y,
               int *ylen, double* z, int *zlen,
               double* clb, int *clblen, double* cub,
               int *cublen, int *kernel, int *dmetric,
               int *n, double* opt_tol, int *opt_max_iters,
               int *ncores, int *gpus, int *ts,
               int *p_grid, int *q_grid, double* globalthetaout);

void mle_exact_non_stat(double* x, int *xlen, double* y,
                        int *ylen, double* z, int *zlen,
                        int *kernel, double* clb, int *clblen,
                        double* cub, int *cublen, int *dmetric,
                        int *n, double* opt_tol, int *opt_max_iters,
                        int *ncores, int *gpus, int *ts,
                        int *p_grid, int *q_grid, double* globalthetaout);

void mle_tlr(double* x, int *xlen, double* y,
             int *ylen, double* z, int *zlen,
             int *kernel, double* clb, int *clblen,
             double* cub, int *cublen, int *tlr_acc,
             int *tlr_maxrank, int *dmetric, int *n,
             double* opt_tol, int *opt_max_iters, int *ncores,
             int *gpus,  int *dts, int *lts, int *p_grid,
             int *q_grid, double* globalthetaout);

void mle_dst(double* x, int *xlen, double* y,
             int *ylen, double* z, int *zlen,
             int *kernel, double* clb, int *clblen,
             double* cub, int *cublen, int *dst_thick,
             int *dmetric, int *n, double* opt_tol,
             int *opt_max_iters, int *ncores, int *gpus,
             int *ts, int *p_grid, int *q_grid,
             double* globalthetaout);

void mle_mp(double* x, int *xlen, double* y,
            int *ylen, double* z, int *zlen,
            int *kernel, double* clb, int *clblen,
            double* cub, int *cublen, int *mp_band,
            int *dmetric, int *n, double* opt_tol,
            int *opt_max_iters, int *ncores, int *gpus,
            int *ts, int *p_grid, int *q_grid,
            double* globalthetaout);

void exact_predict(double* xobs, int *xobs_len, double* yobs,
                   int *yobs_len, double* zobs, int *zobs_len,
                   double* xmiss, int *xmiss_len, double* ymiss,
                   int *ymiss_len, int *nZobs, int *nZmiss,
                   int *kernel, double* est_theta, int *thetalen,
                   int *computation, int *dmetric, int *ncores,
                   int *gpus, int *ts, int *p_grid,
                   int *q_grid, double* globalzpredict);

void exact_mloe_mmom(double* xobs, int *xobs_len, double* yobs,
                     int *yobs_len, double* zobs, int *zobs_len,
                     double* xmiss, int *xmiss_len, double* ymiss,
                     int *ymiss_len, int *nZobs, int *nZmiss,
                     int *kernel, double* est_theta, double* true_theta,
                     int *thetalen, int *dmetric, int *ncores,
                     int *gpus, int *ts, int *p_grid,
                     int *q_grid, double* global_mloe_mmom);

void fisher_general(double* x, double* y, int *n,
                    double* theta, int *thetalen,
                    int *ts, int *dmetric, int *p_grid, int *q_grid,
                    double* globalthetaout);

void  rexageostat_init(int *ncores, int *gpus, int *dts, int *lts);

void rexageostat_finalize();

#endif