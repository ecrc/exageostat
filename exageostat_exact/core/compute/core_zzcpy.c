/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dzcpy.c
 *
 * Copy contents of descriptor to vector.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/exageostatcore.h"



/*******************************************************************************
 *
 *  core_dzcpy - Copy vector r to vector Z (one tile) (double precision).
 *  The routine makes only one pass through the tile Z.
 *  One tile operation.
 *  One vector tile.
 *******************************************************************************
 *
 * @param[out] Z
 *           The m-by-1 matrix on which to store r.
 *
 * @param[in] m
 *          The number of rows in the tile A.
 *
 * @param[in] m0
 *          global row index of the tile A.
 *
 * @param[in] r
 *           double vector.
 *
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_dzcpy(double *Z, int m, int m0, double *r) {

        memcpy(Z, &r[m0], m * sizeof(double));
}




/*******************************************************************************
 *
 *  core_szcpy - Copy vector r to vector Z (one tile) (single precision).
 *  The routine makes only one pass through the tile Z.
 *  One tile operation.
 *  One vector tile.
 *******************************************************************************
 *
 * @param[out] Z
 *           The m-by-1 matrix on which to store r.
 *
 * @param[in] m
 *          The number of rows in the tile A.
 *
 * @param[in] m0
 *          global row index of the tile A.
 *
 * @param[in] r
 *           double vector.
 *
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_szcpy(float *Z,  int m, int m0, float *r) {

        memcpy(Z, &r[m0], m * sizeof(float));
}
