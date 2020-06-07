/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dcmg_nuggets.c
 *
 * Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 1.1.0
 * @author Mary Lai Salvana 
 * @author Sameh Abdulah
 * @date 2020-05-05
 *
 **/
#include "../include/exageostatcore.h"

// This function converts decimal degrees to radians
static double deg2rad(double deg) {
    return (deg * PI / 180);
}
//  This function converts radians to decimal degrees
static double rad2deg(double rad) {
    return (rad * 180 / PI);
}
/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
static double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r)/2);
    v = sin((lon2r - lon1r)/2);
    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

static double calculateDistance(double x1, double y1, double x2, double y2, int distance_metric) {

    if(distance_metric == 1)
        return distanceEarth(x1, y1, x2, y2);
    return  sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}

static double calculateDistance3d(double x1, double y1, double z1,
        double x2, double y2, double z2, int distance_metric) {

    if(distance_metric == 1)
    {
        printf("Great Circle (GC) distance is only valid for 2d\n");
        exit(0);
    }
    return  sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) +  pow((z2 - z1), 2));
}
/***************************************************************************//**
 *
 *  core_dcmg_nuggets - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Matern Kernel).
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
void core_dcmg_nuggets (double *A, int m, int n, int m0, int n0, location  *l1, location *l2, double *localtheta, int distance_metric) {

    //	printf("%f, %f, %f, %f\n", localtheta[0], localtheta[1], localtheta[2], localtheta[3]);
    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0, z0;
    double expr = 0.0;
    double con = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
    con = 1.0/con;
    con = sigma_square * con;
    if(l1->z == NULL || l2->z == NULL )
    {
        for (i = 0; i < m; i++) {
            j0 = n0;
            x0 = l1->x[i0];
            y0 = l1->y[i0];
            for (j = 0; j < n; j++) {
                expr = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric)/localtheta[1];
                if(expr == 0)
                    A[i + j * m] = sigma_square + localtheta[3];
                else
                    A[i + j * m] = con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
                j0++;
            }
            i0++;
        }
    }

    else
    {
        for (i = 0; i < m; i++) {
            j0 = n0;
            x0 = l1->x[i0];
            y0 = l1->y[i0];
            z0 = l1->z[i0];
            for (j = 0; j < n; j++) {
                expr  = calculateDistance3d(x0, y0, z0,
                        l2->x[j0], l2->y[j0], l2->z[j0],  distance_metric);
                if(expr == 0)
                    A[i + j * m] = sigma_square + localtheta[3];
                else
                    A[i + j * m] = con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
                j0++;
            }
            i0++;
        }



    }
}


/***************************************************************************//**
 *
 *  core_scmg_nuggets - Calculate covariance matrix A - single precision.
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
void core_scmg_nuggets (float *A, int m, int n, int m0, int n0, location  *l1, location *l2, double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    float x0, y0;
    float expr = 0.0;
    float con = 0.0;
    float sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
    con = 1.0/con;
    con = sigma_square * con;

    for (i = 0; i < m; i++) {
        j0 = n0;
        x0 = l1->x[i0];
        y0 = l1->y[i0];
        for (j = 0; j < n; j++) {
            expr = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric)/localtheta[1];
            if(expr == 0)
                A[i + j * m] = sigma_square + localtheta[3];
            else
                A[i + j * m] = con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
            j0++;
        }
        i0++;
    }
}







void core_sdcmg_nuggets (double *A, int m, int n, int m0, int n0, location  *l1, location *l2, double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0;
    double expr = 0.0;
    double con = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
    con = 1.0/con;
    con = sigma_square * con;

    for (i = 0; i < m; i++) {
        j0 = n0;
        x0 = l1->x[i0];
        y0 = l1->y[i0];
        for (j = 0; j < n; j++) {
            expr = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric)/localtheta[1];
            if(expr == 0)
                A[i + j * m] = (float)(sigma_square +localtheta[3]);
            else
                A[i + j * m] = (float)(con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr)); // Matern Function
            j0++;
        }
        i0++;
    }
}







/***************************************************************************//**
 *
 *  core_dcmg_nuggets_pow_exp - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Power exp).
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
void core_dcmg_nuggets_pow_exp (double *A, int m, int n,
        int m0, int n0,
        location  *l1, location *l2, 
        double *localtheta, int distance_metric) {

    //Matern kernel
    //localtheta[0] --> variance (sigma),
    //localtheta[1] --> range(beta)
    //localtheta[2] --> smoothness (nu)

    //Power exp kernel
    //localtheta[0] --> variance (sigma),
    //localtheta[1] --> range(beta)
    //localtheta[2] --> range(delta)   0 < delta< 2
    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0;
    double expr  = 0.0;
    double expr1 = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];


    for (i = 0; i < m; i++) {
        j0 = n0;
        x0 = l1->x[i0];
        y0 = l1->y[i0];
        for (j = 0; j < n; j++) {
            expr  = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric);
            expr1 = pow(expr, localtheta[2]);


            if(expr == 0)
                A[i + j * m] = sigma_square + localtheta[3];
            else
                A[i + j * m] = sigma_square *  exp(-(expr1/localtheta[1])) ; // power-exp kernel
            j0++;
        }
        i0++;
    }
}





/***************************************************************************//**
 *
 *  core_scmg_nuggets_pow_exp - Calculate covariance matrix A - single precision.
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
void core_scmg_nuggets_pow_exp (float *A, int m, int n,
        int m0, int n0,
        location  *l1,    location *l2,
        double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    float x0, y0;
    float expr = 0.0;
    float expr1 = 0.0;
    float sigma_square = localtheta[0];// * localtheta[0];


    for (i = 0; i < m; i++) {
        j0 = n0;
        x0 = l1->x[i0];
        y0 = l1->y[i0];
        for (j = 0; j < n; j++) {
            expr  = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric);
            expr1 = pow(expr, localtheta[2]);
            if(expr == 0)
                A[i + j * m] = sigma_square + localtheta[3];
            else
                //A[i + j * m] = con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
                A[i + j * m] = sigma_square *  exp(-(expr1/localtheta[1])) ; // power-exp kernel
            j0++;
        }
        i0++;
    }
}




/***************************************************************************//**
 *
 *  core_sdcmg_nuggets_pow_exp - Calculate covariance matrix A - single precision.
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

void core_sdcmg_nuggets_pow_exp (double *A, int m, int n,
        int m0, int n0, 
        location  *l1, location *l2, 
        double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0;
    double expr = 0.0;
    double expr1 = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];


    for (i = 0; i < m; i++) {
        j0 = n0;
        x0 = l1->x[i0];
        y0 = l1->y[i0];
        for (j = 0; j < n; j++) {
            expr  = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric);
            expr1 = pow(expr, localtheta[2]);
            if(expr == 0)
                A[i + j * m] = (float)(sigma_square + localtheta[3]);
            else
                A[i + j * m] = (float)(sigma_square *  exp(-(expr1/localtheta[1]))); // power-exp kernel
            j0++;
        }
        i0++;
    }
}

