/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_ng_transform.c
 *
 * Transform the tgh measurements vector to Gaussian measurements vectirs, for given parameters. This is required to compute the tgh log-likelihood for that given parameters.
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
 *  core_ng_transform - Transform the tgh measurements vector to Gaussian measurements vectirs, for given parameters.
 *  The routine makes only one pass through the tile A.
 *  One tile operation.
 *******************************************************************************
 *
 * @param[out] Z
 *           The 1-by-n  measurments vector.
 *
 * @param[in] n
 *          The number of cols in the tile A.
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
/***************************************************************************//**
 * 
 *  f - Returns how big is a tgh observation is to a Gaussian observation, for a transformation parameter vector
 *******************************************************************************
 *
 * @param[in] z_non
 *          tgh measurement of a location
 * @param[in] z
 *          possible Gaussian measurement of that location
 * @param[in] xi, omega, g, h
 *          Parameters for the Tukey-gh transformation
 *          xi is location parameter
 *          omega is scale parameter
 *          g is skewness parameter
 *          h is kurtosis parameter
 * returns (z_non - tgh(z,xi,omega,g,h))
 ********************************************************************************
 */
/***************************************************************************//**
 *
 * df - Returns the differentiation of f with respect to z
 *******************************************************************************
 *
 * @param[in] z_non
 *          tgh measurement of a location
 * @param[in] z
 *          possible Gaussian measurement of that location
 * @param[in] xi, omega, g, h
 *          Parameters for the Tukey-gh transformation
 *          xi is location parameter
 *          omega is scale parameter
 *          g is skewness parameter
 *          h is kurtosis parameter
 * returns derivative of (z_non - tgh(z,xi,omega,g,h)) with repect to z
 ********************************************************************************
 */
/***************************************************************************//**
 *
 * newton_raphson - Returns the inverse of the tgh transoformation of a tgh measurement and for given parameters for the Tukey-gh transformation using Newton Raphson algorithm
 *******************************************************************************
 *
 * @param[in] z
 *          tgh measurement of a location
 * @param[in] xi, omega, g, h
 *          Parameters for the Tukey-gh transformation
 *          xi is location parameter
 *          omega is scale parameter
 *          g is skewness parameter
 *          h is kurtosis parameter
 * @param[in] eps
 * 	    Parameter which specifies the maximum error of the solution
 * returns The inverse of the tgh transoformation of a tgh measurement and for given parameters for the Tukey-gh transformation
 ********************************************************************************
 */

static double f(double z_non, double z, double xi, double omega, double g, double h) {
    if (g == 0)
        return z_non - xi - omega * z * exp(0.5 * h * z * z);
    else

        return z_non - xi - (omega * (exp(g * z) - 1) * (exp(0.5 * h * z * z)) / g);
}

static double df(double z, double xi, double omega, double g, double h) {
    if (g == 0)
        return -omega * exp((h * z * z) / 2.0) - omega * h * z * z * exp((h * z * z) / 2.0);
    else
        return -omega * exp(g * z) * exp((h * z * z) / 2.0) -
               (h * z * exp((h * z * z) / 2.0) * (omega * exp(g * z) - omega)) / g;
}

double newton_raphson(double z, double xi, double omega, double g, double h, double eps) {
    int itr, maxmitr;
    double x0, x1, allerr;
    x0 = 0;
    double diff;
    allerr = eps;
    maxmitr = 1000;
    for (itr = 1; itr <= maxmitr; itr++) {
        diff = f(z, x0, xi, omega, g, h) / df(x0, xi, omega, g, h);
        x1 = x0 - diff;
        if (fabs(diff) < allerr)
            return x1;
        x0 = x1;
    }

    return x1;
}

void core_ng_transform(double* Z, double* nan_flag, double* localtheta, int m) {

    double xi = localtheta[2];
    double omega = localtheta[3];
    double g = localtheta[4];
    double h = localtheta[5];

    double eps = 1.0e-5;
    for (int i = 0; i < m; i++)
        Z[i] = newton_raphson(Z[i], xi, omega, g, h, eps);
}
