/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_g_to_ng.c
 *
 * Convert Gaussian random field to Tukey g-and-h random field.
 *
 * @version 1.2.0
 *
 * @author Sagnik Mondal
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/exageostatcore.h"


/***************************************************************************//**
 *
 *  core_g_to_ng - convert Z tile from Gaussian to non-Gaussian.
 *  The routine makes only one pass through the tile A.
 *  One tile operation.    
 *******************************************************************************
 *
 * @param[out] Z
 *           The 1-by-n vector on which to compute the covariance matrix.
 *
 * @param[in] n
 *          The number of cols in the tile Z. 
 *
 * @param[in] localtheta[0],localtheta[1]
 *          Parameter vector that is used to generate the output covariance matrix.
 *
 * @param[in] localtheta[2],localtheta[3],localtheta[4],localtheta[5]
 *          Parameter vector that is used in Tukey-gh transformation
 *          localtheta[2] is location parameter
 *          localtheta[3] is scale parameter
 *          localtheta[4] is skewness parameter
 *          localtheta[5] is kurtosis parameter
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_g_to_ng(double* Z, double* localtheta, int m) {

    double xi = localtheta[2];
    double omega = localtheta[3];
    double g = localtheta[4];
    double h = localtheta[5];

    int i;
    if (h < 0) {
        printf("The kurtosis parameter cannot be negative");
        return;
    }
    if (g == 0) {
        for (i = 0; i < m; i++)
            Z[i] = xi + omega * Z[i] * (exp(0.5 * h * pow(Z[i], 2)));
    } else {
        for (i = 0; i < m; i++)
            Z[i] = xi + omega * (exp(g * Z[i]) - 1) * (exp(0.5 * h * pow(Z[i], 2))) / g;
    }
}