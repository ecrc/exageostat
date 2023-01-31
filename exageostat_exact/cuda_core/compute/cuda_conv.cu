/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file cudaconv.c
 *
 * Cuda datatypes conversion.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/

#define CHUNKSIZE 32

#include <cublas.h>
#include <stdio.h>
#include "../include/exageostatcudacore.h"

__global__ void
float2double_array_kernel(int nrows, int ncols, const float *F, int ldf, double* H, int ldh, cublasOperation_t transa) {
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int idx = blockIdx.x * blockDim.x + tx;
    const int idy = blockIdx.y * blockDim.y + ty;

    if (idx >= nrows || idy >= ncols) {
        return;
    }
    if (transa == CUBLAS_OP_N)
        H[idy * ldh + idx] = (double) F[idy * ldf + idx];
    else
        H[idx * ldh + idy] = (double) F[idy * ldf + idx];
}

void float2double_array(int nrows, int ncols, const float *F, int ldf, double* H, int ldh, cublasOperation_t transa,
                        cudaStream_t stream) {
    int nBlockx = (nrows + CHUNKSIZE - 1) / CHUNKSIZE;
    int nBlocky = (ncols + CHUNKSIZE - 1) / CHUNKSIZE;
    dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
    dim3 dimGrid(nBlockx, nBlocky);
    float2double_array_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, F, ldf, H, ldh, transa);
}

__global__ void
double2float_array_kernel(int nrows, int ncols, const double* H, int ldh, float *F, int ldf, cublasOperation_t transa) {
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int idx = blockIdx.x * blockDim.x + tx;
    const int idy = blockIdx.y * blockDim.y + ty;

    if (idx >= nrows || idy >= ncols) {
        return;
    }

    if (transa == CUBLAS_OP_N)
        F[idy * ldf + idx] = __double2float_rn(
                H[idy * ldh + idx]); //Convert a double to a float in round-to-nearest-even mode.
    else
        F[idx * ldf + idy] = __double2float_rn(
                H[idy * ldh + idx]); //Convert a double to a float in round-to-nearest-even mode
}

void double2float_array(int nrows, int ncols, const double* H, int ldh, float *F, int ldf, cublasOperation_t transa,
                        cudaStream_t stream) {

    int nBlockx = (nrows + CHUNKSIZE - 1) / CHUNKSIZE;
    int nBlocky = (ncols + CHUNKSIZE - 1) / CHUNKSIZE;
    dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
    dim3 dimGrid(nBlockx, nBlocky);
    double2float_array_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, H, ldh, F, ldf, transa);
}