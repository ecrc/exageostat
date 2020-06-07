/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dmdet.c
 *
 * Calculate determinant of a given triangular matrix (A).
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/

#include "../include/exageostatcore.h"


/***************************************************************************//**
 *
 *  core_dmdet - Calculate the determinant of the matrix A (double precision).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.
 *******************************************************************************
 *
 * @param[out] A
 *           The m-by-n matrix on which to calculate the determinant.
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
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
double core_dmdet (double *A, int m, int n, int m0, int n0) {

	int i;
	double res = 0.0;
	for (i = 0; i < m; i++) {
		if(A[i + i * m] > 0)
			res += log(A[i + i * m]);
	//	printf("(in loop) %s, %.7g,  %.7g \n", __func__, log(A[i + i * m]),  A[i + i * m]);
	}

	//printf("%s, %f\n", __func__, res); 

	return res;
}





/***************************************************************************//**
 *
 *  core_smdet - Calculate the determinant of the matrix A (single precision).
 *  The routine makes only one pass through the tile A.
 *  One tile operation.
 *******************************************************************************
 *
 * @param[out] A
 *           The m-by-n matrix on which to calculate the determinant.
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
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
float core_smdet (float * A, int m, int n, int m0, int n0) {

	int i;
	float res = 0.0;
	for (i = 0; i < m; i++) {
		if(A[i + i * m] > 0)
			res += log(A[i + i * m]);
	}

	return res;
}




