/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dcmg.c
 *
 * Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 1.1.0
 * @author Mary Lai Salvana 
 * @author Sameh Abdulah
 * @date 2020-03-29
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

static double calculateDistance( location* l1, location* l2, int l1_index,
        int l2_index, int distance_metric, int z_flag) {

    double z1, z2;
    double x1=l1->x[l1_index];
    double y1=l1->y[l1_index];
    double x2=l2->x[l2_index];
    double y2=l2->y[l2_index];
    if(l1->z == NULL || l2->z == NULL || z_flag == 1)
    {
        if(distance_metric == 1)
            return distanceEarth(x1, y1, x2, y2);
        return  sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
    }
    else
    {
        if(distance_metric == 1)
        {
            printf("Great Circle (GC) distance is only valid for 2d\n");
            exit(0);
        }
        z1 = l1->z[l1_index];
        z2 = l2->z[l2_index];
        return  sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) +  pow((z2 - z1), 2));
    }	
}

/*
   static double calculateDistance3d(double x1, double y1, double z1, 
   double x2, double y2, double z2, int distance_metric) {

   if(distance_metric == 1)
   {
   printf("Great Circle (GC) distance is only valid for 2d\n");
   exit(0);
   }
   return  sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) +  pow((z2 - z1), 2));
   }
   */
/***************************************************************************//**
 *
 *  core_dcmg - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Matern Kernel).
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
void core_dcmg (double *A, int m, int n,
        int m0, int n0, location  *l1,
        location *l2, double *localtheta, int distance_metric) {

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

    for (i = 0; i < m; i++) {
        j0 = n0;
        for (j = 0; j < n; j++) {
            expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0)/localtheta[1];
            if(expr == 0)
                A[i + j * m] = sigma_square /*+ 1e-4*/;
            else
                A[i + j * m] = con*pow(expr, localtheta[2])
                    * gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
            j0++;
        }
        i0++;
    }

}


/***************************************************************************//**
 *
 *  core_scmg - Calculate covariance matrix A - single precision.
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
void core_scmg (float *A, int m, int n,
        int m0, int n0, location  *l1,
        location *l2, double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    float x0, y0, z0;
    float expr = 0.0;
    float con = 0.0;
    float sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
    con = 1.0/con;
    con = sigma_square * con;


    for (i = 0; i < m; i++) {
        j0 = n0;
        for (j = 0; j < n; j++) {
            expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
            if(expr == 0)
                A[i + j * m] = sigma_square /*+ 1e-4*/;
            else
                A[i + j * m] = con*pow(expr, localtheta[2])
                    * gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
            j0++;
        }
        i0++;
    }
}







void core_sdcmg (double *A, int m, int n, 
        int m0, int n0, location  *l1,
        location *l2, double *localtheta, int distance_metric) {

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


    for (i = 0; i < m; i++) {
        j0 = n0;
        for (j = 0; j < n; j++) {
            expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0)/localtheta[1];
            if(expr == 0)
                A[i + j * m] = (float)(sigma_square /*+ 1e-4*/);
            else
                A[i + j * m] = (float)(con*pow(expr, localtheta[2])
                        * gsl_sf_bessel_Knu(localtheta[2],expr)); // Matern Function
            j0++;
        }
        i0++;
    }

}




/***************************************************************************//**
 *
 *  core_dcmg_spacetime_matern - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Matern Kernel).
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
void core_dcmg_spacetime_matern (double *A, int m, int n,
        int m0, int n0, location  *l1,
        location *l2, double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0, z0, z1;
    double expr = 0.0, expr1 = 0.0, expr2 = 0.0, expr3 = 0.0, expr4 = 0.0;
    double con = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
    con = 1.0/con;
    con = sigma_square * con;

    for (i = 0; i < m; i++) {
        j0 = n0;
        z0 = l1->z[i0];
        for (j = 0; j < n; j++) {
            z1 = l2->z[j0];
            expr = calculateDistance(l1, l2, i0, j0, distance_metric, 1)/localtheta[1];
            expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * localtheta[4])/localtheta[3] + 1, localtheta[5] / 2);
            expr3 = expr / expr2;
            expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * localtheta[4]) / localtheta[3] + 1, localtheta[5] + localtheta[6]);

            if(expr == 0)
                A[i + j * m] = sigma_square / expr4 /*+ 1e-4*/;
            else
                A[i + j * m] = con*pow(expr3, localtheta[2])
                    * gsl_sf_bessel_Knu(localtheta[2],expr3) / expr4; // Matern Function
            //printf("%f,%f,%f,%f,%f\n", calculateDistance(l1, l2, i0, j0, distance_metric, 1), z0, z1, sqrt(pow(z0 - z1, 2)), A[i + j * m]);
            j0++;
        }
        i0++;
    }

}


/***************************************************************************//**
 *
 *  core_dcmg_bivariate_flexible - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Flexible Bivariate Matern Kernel).
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
 *          localtheta[0] is the variance of Z1.
 *          localtheta[1] is the variance of Z2.
 *          localtheta[2] is the scale of Z1.
 *          localtheta[3] is the scale of Z2.
 *          localtheta[4] is the first latent scale parameter for the cross of Z1 and Z2.
 *          localtheta[5] is the second latent scale parameter for the cross of Z1 and Z2.
 *          localtheta[6] is the smoothness of Z1.
 *          localtheta[7] is the smoothness of Z2.
 *          localtheta[8] is the first latent smoothness parameter for the cross of Z1 and Z2.
 *          localtheta[9] is the second latent smoothness parameter for the cross of Z1 and Z2.
 *          localtheta[10] is the first latent variance parameter for the cross of Z1 and Z2.
 *          localtheta[11] is the second latent variance parameter for the cross of Z1 and Z2.
 *          localtheta[12] is the third latent variance parameter for the cross of Z1 and Z2.
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/

void core_dcmg_bivariate_flexible (double *A, int m, int n,
        int m0, int n0, location  *l1,
        location *l2, double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0;
    double expr1 = 0.0, expr2 = 0.0, expr12 = 0.0;
    double con1 = 0.0, con2 = 0.0, con12 = 0.0, scale12 = 0.0, rho = 0.0, nu12 = 0.0;
    double scale1 = localtheta[0], scale2 = localtheta[1], nu1 = localtheta[4], nu2 = localtheta[5];

    scale12 =  pow(0.5 * (pow(scale1, -2) + pow(scale2, -2)) + localtheta[2] * (1 - localtheta[3]), -0.5) ;

    nu12 = 0.5 * (nu1 + nu2) + localtheta[6] * (1 - localtheta[7]) ;

    rho = localtheta[8] * localtheta[9] * localtheta[10] * 
        pow(scale12,  2 * localtheta[6] + (nu1 + nu2))
        * tgamma(0.5 * (nu1 + nu2) + 1) * tgamma(nu12) / tgamma(nu12 + 1);

    localtheta[0] = localtheta[8] * localtheta[8] * 
        pow(scale1, 2*localtheta[6] + nu1 + nu1) * tgamma(nu1);

    localtheta[1] = localtheta[9] * localtheta[9] *  
        pow(scale2, 2*localtheta[6] + nu2 + nu2) * tgamma(nu2);

    con1 = pow(2,(nu1-1)) * tgamma(nu1);
    con1 = 1.0/con1;
    con1 = localtheta[0] * con1;

    con2 = pow(2,(nu2-1)) * tgamma(nu2);
    con2 = 1.0/con2;
    con2 = localtheta[1] * con2;    

    con12 = pow(2,(nu12-1)) * tgamma(nu12);
    con12 = 1.0/con12;
    con12 = rho * con12;    

    i0/=2;
    for (i = 0; i < m; i+=2) {
        j0 = n0/2;
        for (j = 0; j < n; j+=2) {
            expr1  = calculateDistance(l1, l2, i0, j0, distance_metric, 0)/scale1;
            expr2  = calculateDistance(l1, l2, i0, j0, distance_metric, 0)/scale2;
            expr12 = calculateDistance(l1, l2, i0, j0, distance_metric, 0)/scale12;

            if(expr1 == 0){
                A[i + j * m] = localtheta[0] ;
                A[(i + 1) + j * m] = A[i + (j + 1) * m] = rho;
                A[(i + 1) + (j + 1) * m] = localtheta[1] ;
            }
            else{
                A[i + j * m] = con1 * pow(expr1, nu1)
                    * gsl_sf_bessel_Knu(nu1, expr1); 
                A[(i + 1) + j * m] = A[i + (j + 1) * m] = con12 
                    * pow(expr12, nu12) * gsl_sf_bessel_Knu(nu12, expr12);  
                A[(i + 1) + (j + 1) * m] = con2 * pow(expr2, nu2) 
                    * gsl_sf_bessel_Knu(nu2, expr2);
            }
            j0++;
        }
        i0++;
    }
}



/***************************************************************************//**
 *
 *  core_dcmg_bivariate_parsimonious - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Parsimonious Bivariate Matern Kernel).
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
 *          localtheta[0] is the variance of Z1.
 *          localtheta[1] is the variance of Z2.
 *          localtheta[2] is the common scale parameter for Z1 and Z2.
 *          localtheta[3] is the smoothness of Z1.
 *          localtheta[4] is the smoothness of Z2.
 *          localtheta[5] is the correlation parameter of Z1 and Z2.
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_dcmg_bivariate_parsimonious (double *A, int m, int n,
        int m0, int n0, location  *l1,
        location *l2, double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0;
    double expr = 0.0;
    double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;
    //	double localtheta[0] = localtheta[0], localtheta[1] = localtheta[1];

    con1 = pow(2,(localtheta[3]-1)) * tgamma(localtheta[3]);
    con1 = 1.0/con1;
    con1 = localtheta[0] * con1;

    con2 = pow(2,(localtheta[4]-1)) * tgamma(localtheta[4]);
    con2 = 1.0/con2;
    con2 = localtheta[1] * con2;

    //The average
    nu12 = 0.5 * (localtheta[3] + localtheta[4]);    

    rho = localtheta[5] * sqrt( (tgamma(localtheta[3] + 1)*tgamma(localtheta[4] + 1)) /
            (tgamma(localtheta[3]) * tgamma(localtheta[4])) ) *
        tgamma(nu12) / tgamma(nu12 + 1);



    con12 = pow(2,(nu12-1)) * tgamma(nu12);
    con12 = 1.0/con12;
    con12 = rho * sqrt(localtheta[0] * localtheta[1]) * con12;    

    i0/=2;
    for (i = 0; i < m-1; i+=2) {
        j0 = n0/2;
        for (j = 0; j < n-1; j+=2) {
            expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0)/localtheta[2];

            if(expr == 0){
                A[i + j * m] = localtheta[0] ;
                A[(i + 1) + j * m] = A[i + (j + 1) * m] = rho 
                    * sqrt(localtheta[0] * localtheta[1]) ;
                A[(i + 1) + (j + 1) * m] = localtheta[1] ;
            }
            else{
                A[i + j * m] = con1 * pow(expr, localtheta[3])
                    * gsl_sf_bessel_Knu(localtheta[3], expr); 
                A[(i + 1) + j * m] = A[i + (j + 1) * m] = con12 * pow(expr, nu12) 
                    * gsl_sf_bessel_Knu(nu12, expr);  
                A[(i + 1) + (j + 1) * m] = con2 * pow(expr, localtheta[4])
                    * gsl_sf_bessel_Knu(localtheta[4], expr);
            }
            // printf ("===%d, %d, %d, %d (%f)\n" , i, j, i0, j0, expr);
            j0++;
        }
        i0++;
    }
}

/*
   void core_dcmg_bivariate_parsimonious2(double *A, int m, int n, int m0, int n0, location  *l1, location *l2, double *localtheta, int distance_metric, int size) {
   int i, j;
   int i0 = m0;
   int j0 = n0;
   double x0, y0;
   double expr = 0.0;
   double con = 0.0; double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;
   double sigma_square = 0.0; double nu = 0.0;
   double localtheta[0]2 = 0.0; 
   double localtheta[0]  = localtheta[0];
   double localtheta[1]  = localtheta[1];
   double beta    = localtheta[2];
   double nu1     = localtheta[3];
   double nu2     = localtheta[4];
   double  a      = localtheta[5];
   con1 = pow(2,(localtheta[3]-1)) * tgamma(localtheta[3]);
   con1 = 1.0/con1;
   con1 = localtheta[0] * con1;
   con2 = pow(2,(localtheta[4]-1)) * tgamma(localtheta[4]);
   con2 = 1.0/con2;
   con2 = localtheta[1] * con2;
//The average
nu12 = 0.5 * (localtheta[3] + localtheta[4]);
rho = localtheta[5] * sqrt( (tgamma(localtheta[3] + 1)*tgamma(localtheta[4] + 1)) /
(tgamma(localtheta[3]) * tgamma(localtheta[4])) ) *
tgamma(nu12) / tgamma(nu12 + 1);
con12 = pow(2,(nu12-1)) * tgamma(nu12);
con12 = 1.0/con12;
con12 = rho * sqrt(localtheta[0] * localtheta[1]) * con12;
if(i0<size/2 && j0<size/2)
{
sigma_square = localtheta[0];
nu = nu1;
con =con1;
}
else  if(i0>=size/2 && j0>=size/2)
{
sigma_square = localtheta[1];
nu = nu2;
con=con2;
}
else
{
nu12 = 0.5 * (nu1 + nu2);
rho = a * sqrt( (tgamma(nu1 + 1)*tgamma(nu2 + 1)) /
(tgamma(nu1) * tgamma(nu2)) ) *
tgamma(nu12) / tgamma(nu12 + 1);
localtheta[0]2 = rho * sqrt(localtheta[0]*localtheta[1]) ;
sigma_square = localtheta[0]2;
nu = nu12;
con =con12;
}
for (i = 0; i < m; i++) {
j0 = n0;
if(j0>=size/2)
j0-=size/2;
if(i0 >=size/2)
i0-=size/2;
x0 = l1->x[i0];
y0 = l1->y[i0];
for (j = 0; j < n; j++) {
expr = calculateDistance(x0, y0, l2->x[j0], l2->y[j0], distance_metric)/beta;
if(expr == 0)
if (i+m0>=size/2.0 && j+n0>=size/2.0)
A[i + j * m] = localtheta[1] ;
else if (i+m0<size/2.0 && j+n0<size/2.0)
A[i + j * m] = localtheta[0] ;
else
A[i + j * m] = localtheta[0]2 ;
else
{
    if (i+m0>=size/2.0 && j+n0>=size/2.0)
        A[i + j * m] = con2*pow(expr, nu2)*gsl_sf_bessel_Knu(nu2,expr); // Matern Function
    else if (i+m0<size/2.0 && j+n0<size/2)
        A[i + j * m] = con1*pow(expr, nu1)*gsl_sf_bessel_Knu(nu1,expr); // Matern Function                
    else 
        A[i + j * m] = con12*pow(expr, nu12)*gsl_sf_bessel_Knu(nu12,expr); // Matern Function
}
j0++;
}
i0++;
}
}
*/

/***************************************************************************//**
 *
 *  core_dcmg_pow_exp - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Power exp).
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
void core_dcmg_pow_exp (double *A, int m, int n,
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
    double x0, y0, z0;
    double expr  = 0.0;
    double expr1 = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    for (i = 0; i < m; i++) {
        j0 = n0;
        for (j = 0; j < n; j++) {
            expr  = calculateDistance(l1, l2, i0, j0, distance_metric, 0);
            expr1 = pow(expr, localtheta[2]);


            if(expr == 0)
                A[i + j * m] = sigma_square /*+ 1e-4*/;
            else
                A[i + j * m] = sigma_square *  exp(-(expr1/localtheta[1])) ; // power-exp kernel
            j0++;
        }
        i0++;
    }
}



/***************************************************************************//**
 *
 *  core_scmg_pow_exp - Calculate covariance matrix A - single precision.
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
void core_scmg_pow_exp (float *A, int m, int n,
        int m0, int n0,
        location  *l1,    location *l2,
        double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    float x0, y0, z0;
    float expr = 0.0;
    float expr1 = 0.0;
    float sigma_square = localtheta[0];// * localtheta[0];


    for (i = 0; i < m; i++) {
        j0 = n0;
        for (j = 0; j < n; j++) {
            expr  = calculateDistance(l1, l2,
                    i0, j0, distance_metric, 0);
            expr1 = pow(expr, localtheta[2]);
            if(expr == 0)
                A[i + j * m] = sigma_square /*+ 1e-4*/;
            else
                A[i + j * m] = sigma_square *  exp(-(expr1/localtheta[1])) ; // power-exp kernel
            j0++;
        }
        i0++;
    }

}


/***************************************************************************//**
 *
 *  core_sdcmg_pow_exp - Calculate covariance matrix A - single precision.
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

void core_sdcmg_pow_exp (double *A, int m, int n,
        int m0, int n0, 
        location  *l1, location *l2, 
        double *localtheta, int distance_metric) {

    int i, j;
    int i0 = m0;
    int j0 = n0;
    double x0, y0, z0;
    double expr = 0.0;
    double expr1 = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    for (i = 0; i < m; i++) {
        j0 = n0;
        for (j = 0; j < n; j++) {
            expr  = calculateDistance(l1, l2, i0, j0, distance_metric, 0);
            expr1 = pow(expr, localtheta[2]);
            if(expr == 0)
                A[i + j * m] = (float)(sigma_square /*+ 1e-4*/);
            else
                A[i + j * m] = (float)(sigma_square *  exp(-(expr1/localtheta[1]))); // power-exp kernel
            j0++;
        }
        i0++;
    }
}


/******************/
//void core_dcmg_deformation (double *A, int m, int n,
//		int m0, int n0, location  *l1,
//		location *l2, double *localtheta, int distance_metric) {

//	int i, j;
//	int i0 = m0;
//	int j0 = n0;
//	double x0, y0, x1, y1, deform_x0, deform_x1, deform_y0, deform_y1;
//	double expr = 0.0;
//	double con = 0.0;
//	double sigma_square = localtheta[0], pt_src1 = localtheta[3], pt_src2 = localtheta[4];// * localtheta[0];

//	con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
//	con = 1.0/con;
//	con = sigma_square * con;

//	for (i = 0; i < m; i++) {
//		j0 = n0;
//		x0 = l1->x[i0];
//		y0 = l1->y[i0];
//		deform_x0 = pt_src1 + (x0 - pt_src1) * calculateDistance(x0, y0, pt_src1, pt_src2, distance_metric);
//		deform_y0 = pt_src2 + (y0 - pt_src2) * calculateDistance(x0, y0, pt_src1, pt_src2, distance_metric);

//		for (j = 0; j < n; j++) {

//			x1 = l1->x[j0];
//			y1 = l1->y[j0];

//			deform_x1 = pt_src1 + (x1 - pt_src1) * calculateDistance(x1, y1, pt_src1, pt_src2, distance_metric);
//			deform_y1 = pt_src2 + (y1 - pt_src2) * calculateDistance(x1, y1, pt_src1, pt_src2, distance_metric);
//			expr = calculateDistance(deform_x0, deform_y0, deform_x1, deform_y1, distance_metric)/localtheta[1];

//			if(expr == 0)
//				A[i + j * m] = sigma_square /*+ 1e-4*/;
//			else
//				A[i + j * m] = con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
//			j0++;
//		}
//		i0++;
//	}
//}
//*/
