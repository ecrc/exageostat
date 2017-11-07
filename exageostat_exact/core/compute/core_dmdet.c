/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/

#include "../include/exageostatcore.h"


/***************************************************************************//**
 *
 *  core_dmdet - Calculate the determinant of the matrix A.
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
double core_dmdet (double * A, int m, int n, int m0, int n0) {

        int i;
        double res = 0.0;
        for (i = 0; i < m; i++) {
                res += log(A[i + i * m]);
        }

        return res;
}
