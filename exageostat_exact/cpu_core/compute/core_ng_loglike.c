/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_ng_loglike.c
 *
 * Calculate the loglikelihood of non-Gaussian MLE.
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
 *  core_ng_loglike - Calculate the loglikelihood of non-Gaussian MLE.
 *  The routine makes only one pass through the tile A.
 *  One tile operation.	
 *******************************************************************************
 *
 * @param[out] Z
 *           Returns the extra terms which are required to be added to the Gaussian log-likelihood
 *           to calculate the tgh log-likelihood
 *
 * @param[in] n
 *          The number of cols in the tile A.
 *
 * @param[in] localtheta[0], localtheta[1]
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
double core_ng_loglike(double* Z, double* localtheta, int m) {

    double xi = localtheta[2];
    double omega = localtheta[3];
    double g = localtheta[4];
    double h = localtheta[5];

    int i;
    double sum = 0;
    if (h < 0) {
        printf("kurtosis parameter cannot be negative");
        exit(1);
    }
    for (i = 0; i < m; i++) {
        if (g == 0)
            sum += log(1 + h * pow(Z[i], 2)) + 0.5 * h * pow(Z[i], 2);
        else {
            sum += log(exp(g * Z[i]) + (exp(g * Z[i]) - 1) * h * Z[i] / g) + 0.5 * h * pow(Z[i], 2);
        }
    }
    return (sum);
}
