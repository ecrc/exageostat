/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/

#define CHUNKSIZE 32

#include <cublas.h>
#include <stdio.h>
#include "../include/exageostatcudacore.h"


__global__ void dcmg_array_kernel(double *A, int m, int n, int m0,
        int n0, double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
        double localtheta0, double localtheta1, double localtheta2, int distance_metric)
{

    const int tx  = threadIdx.x;
    const int ty  = threadIdx.y;
    const int idx = blockIdx.x * blockDim.x + tx;
    const int idy = blockIdx.y * blockDim.y + ty;

    if(idx>=m || idy >=n){return;}

    //double x0, y0;
    double expr  = 0.0;
    double expr1 = 0.0;

    double sigma_square = localtheta0;// * localtheta[0];

    expr = sqrt(pow((l2_x_cuda[idx] - l1_x_cuda[idy]), 2) +
            pow((l2_y_cuda[idx] - l1_y_cuda[idy]), 2));

    expr1 = pow(expr, localtheta2);
    if(expr == 0)
        A[idx + idy * m] = sigma_square /*+ 1e-4*/;
    else
        A[idx + idy * m] = sigma_square *  exp(-(expr1/localtheta1)); // power-exp kernel



}

void dcmg_array( double *A, int m, int n, int m0,
        int n0, double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
        double *localtheta, int distance_metric, cudaStream_t stream){

    int nBlockx= (m+CHUNKSIZE-1)/CHUNKSIZE;
    int nBlocky= (n+CHUNKSIZE-1)/CHUNKSIZE;
    dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
    dim3 dimGrid(nBlockx,nBlocky);


    dcmg_array_kernel<<<dimGrid,dimBlock,0,stream>>>(A, m, n, m0, n0, l1_x_cuda, l1_y_cuda, l2_x_cuda, l2_y_cuda, localtheta[0],localtheta[1],localtheta[2], distance_metric);



}
