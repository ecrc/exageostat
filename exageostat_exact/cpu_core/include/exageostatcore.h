/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file exageostatcore.h
 *
 * Core functions header file.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#ifndef _EXAGEOSTATCORE_H_
#define _EXAGEOSTATCORE_H_

#include "../../../misc/include/MLE_misc.h"
#include "../../../misc/include/flat_file.h"

//Generate the covariance matrix.
void core_scmg(float *A, int m, int n,
               int m0, int n0,
               location *l1, location *l2,
               double* localtheta, int distance_metric);

void core_dcmg(double* A, int m, int n,
               int m0, int n0,
               location *l1, location *l2,
               double* localtheta, int distance_metric);


void core_sdcmg(double* A, int m, int n,
                int m0, int n0,
                location *l1, location *l2,
                double* localtheta, int distance_metric);


void core_scmg_pow_exp(float *A, int m, int n,
                       int m0, int n0,
                       location *l1, location *l2,
                       double* localtheta, int distance_metric);

void core_dcmg_pow_exp(double* A, int m, int n,
                       int m0, int n0,
                       location *l1, location *l2,
                       double* localtheta, int distance_metric);


void core_sdcmg_pow_exp(double* A, int m, int n,
                        int m0, int n0,
                        location *l1, location *l2,
                        double* localtheta, int distance_metric);

void core_dcmg_bivariate_parsimonious(double* A, int m, int n,
                                      int m0, int n0, location *l1,
                                      location *l2, double* localtheta, int distance_metric);

void core_dcmg_bivariate_parsimonious2(double* A, int m, int n,
                                       int m0, int n0, location *l1,
                                       location *l2, double* localtheta, int distance_metric, int size);

void core_dcmg_bivariate_flexible(double* A, int m, int n,
                                  int m0, int n0, location *l1,
                                  location *l2, double* localtheta, int distance_metric);

float core_smdet(float *A, int m, int n,
                 int m0, int n0);

double core_dmdet(double* A, int m, int n,
                  int m0, int n0);

void core_szcpy(float *Z, int m,
                int m0, float *r);

void core_dzcpy(double* Z, int m,
                int m0, double* r);

float core_sdotp(float *Z, float *dotproduct,
                 int n);

double core_ddotp(double* Z, double* dotproduct,
                  int n);

void core_dlag2s(int m, int n,
                 const double* A, int lda,
                 float *B, int ldb);

void core_slag2d(int m, int n,
                 const float *A, int lda,
                 double* B, int ldb);

void core_sprint(float *A,
                 int m, int n,
                 int m0, int n0);

void core_dprint(double* A,
                 int m, int n,
                 int m0, int n0);

void core_dcmg_nono_stat(double* A, int m, int n,
                         int m0, int n0, location *l1,
                         location *l2, location *lm, double* localtheta,
                         int distance_metric);

void core_dcmg_matern_dsigma_square(double* A, int m, int n,
                                    int m0, int n0, location *l1,
                                    location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_dnu(double* A, int m, int n,
                          int m0, int n0, location *l1,
                          location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_dbeta(double* A, int m, int n,
                            int m0, int n0, location *l1,
                            location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_ddsigma_square(double* A, int m, int n);

void core_dcmg_matern_ddsigma_square_beta(double* A, int m, int n,
                                          int m0, int n0, location *l1,
                                          location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_ddsigma_square_nu(double* A, int m, int n,
                                        int m0, int n0, location *l1,
                                        location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_ddbeta_beta(double* A, int m, int n,
                                  int m0, int n0, location *l1,
                                  location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_ddbeta_nu(double* A, int m, int n,
                                int m0, int n0, location *l1,
                                location *l2, double* localtheta, int distance_metric);

double core_dtrace(double* A, int m, int n,
                   int m0, int n0, double* trace);

double core_ng_loglike(double* Z, double* localtheta,
                       int m);

void core_ng_transform(double* Z, double* nan_flag, double* localtheta,
                       int m);

void core_g_to_ng(double* Z, double* localtheta,
                  int m);

double core_dtrace(double* A, int m, int n,
                   int m0, int n0, double* trace);

void core_dcmg_nuggets(double* A, int m, int n,
                       int m0, int n0, location *l1,
                       location *l2, double* localtheta, int distance_metric);

void core_dcmg_spacetime_bivariate_parsimonious(double* A, int m, int n,
                                                int m0, int n0, location *l1,
                                                location *l2, double* localtheta, int distance_metric);

void core_dcmg_non_stat(double* A, int m, int n, int m0,
                        int n0, location *l1, location *l2, double* localtheta, int distance_metric);

void core_dcmg_spacetime_matern(double* A, int m, int n,
                                int m0, int n0, location *l1,
                                location *l2, double* localtheta, int distance_metric);

void core_dcmg_matern_ddnu_nu(double* A, int m, int n, int m0, int n0, location *l1, location *l2, double* localtheta,
                              int distance_metric);

void core_ng_dcmg(double* A, int m, int n,
                  int m0, int n0, location *l1,
                  location *l2, double* localtheta, int distance_metric);

void core_ng_exp_dcmg(double* A, int m, int n,
                      int m0, int n0, location *l1,
                      location *l2, double* localtheta, int distance_metric);


void core_dcmg_trivariate_parsimonious(double* A, int m, int n,
                                       int m0, int n0, location *l1,
                                       location *l2, double* localtheta, int distance_metric);

#endif