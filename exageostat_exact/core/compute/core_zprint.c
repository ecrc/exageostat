/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dprint.c
 *
 * Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-26
 *
 **/
#include "../include/exageostatcore.h"

/***************************************************************************//**
 *
 *  core_dprint - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Matern Kernel).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.	
 *******************************************************************************
 *
 * @param[out] A
 *           The m-by-n matrix on which to compute the covariance matrix.
 *
 * @param[in] m
 *          The number of rows in the tile A. 
 *
 * @param[in] n
 *          The number of cols in the tile A. 
 *
 * @param[in] m0
 *          Global row index of the tile A.
 *
 * @param[in] n0
 *          Global col index of the tile A.
 *
 * @param[in] l1
 *          Location struct of the first input.
 *
 * @param[in] l2
 *          Location struct of the second input.
 *
 * @param[in] localtheta
 *          Parameter vector that is used to generate the output covariance matrix.
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_dprint (double *A, int m, int n, int m0, int n0) {

	int i, j;

	if (n0 == 0)
		printf("\n");

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			printf("%f - ", A[i + j * m]);
		}
	}

}


/***************************************************************************//**
 *
 *  core_sprint - Calculate covariance matrix A - single precision.
 *  The routine makes only one pass through the tile A.
 *  One tile operation.
 *******************************************************************************
 *
 * @param[out] A
 *           The m-by-n matrix on which to compute the covariance matrix.
 *
 * @param[in] m
 *          The number of rows in the tile A.
 *
 * @param[in] n
 *          The number of cols in the tile A.
 *
 * @param[in] m0
 *          global row index of the tile A.
 *
 * @param[in] n0
 *           global col index of the tile A.
 *
 * @param[in] l1
 *          location struct of the first input.
 *
 * @param[in] l2
 *          location struct of the second input.
 *
 * @param[in] localtheta
 *           parameter vector that should be used to generate the output covariance matrix
 *
 * @param[in] distance_metric
 *           distance metric "euclidean Distance (ED"" or "Great Circle Distance (GCD)"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_sprint (float *A, int m, int n, int m0, int n0) {

	int i, j;

	if (n0 == 0)
		printf("\n");

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			printf("%f - ", A[i + j * m]);

		}
	}

}

