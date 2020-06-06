/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file exageostatcudacore.h
 *
 * CUDA core functions header file.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#ifndef _EXAGEOSTATCUDACORE_H_
#define _EXAGEOSTATCUDACORE_H_
#include "cudablas.h"
#include "lapacke.h"
#include "../../misc/include/MLE_misc.h"
#include "../../misc/include/flat_file.h"
int cuda_dsconv(int m, int n,
        cuDoubleComplex *A, int lda,
        cuFloatComplex *B, int ldb, cublasHandle_t handle);


int cuda_sdconv(int m, int n,
        cuFloatComplex *A, int lda,
        cuDoubleComplex *B, int ldb, cublasHandle_t handle);


void cuda_dcmg( double *A, int m, int n, int m0,
        int n0, double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
        double *localtheta, int distance_metric, cublasHandle_t handle);
#ifdef __cplusplus
extern "C"
{
#endif
    void double2float_array(int nrows, int ncols,
            const double *H, int ldh,
            float *F, int ldf,
            cublasOperation_t transa, cudaStream_t stream);

    void float2double_array(int nrows, int ncols,
            const float *F, int ldf,
            double *H, int ldh ,
            cublasOperation_t transa, cudaStream_t stream);


    void dcmg_array(double *A, int m, int n,
            int m0, int n0,
            double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
            double *localtheta, int distance_metric,
            cudaStream_t stream);

#ifdef __cplusplus
}
#endif
#endif


