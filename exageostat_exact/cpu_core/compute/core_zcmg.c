/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
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
 * @version 1.2.0
 * @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/exageostatcore.h"
#include "assert.h"

#define pow_e (2.71828182845904)

/**
 * This function converts decimal degrees to radians
 * @param deg decimal degree
 */
static double deg2rad(double deg) {
	return (deg * PI / 180);
}

/**
 * This function converts radians to decimal degrees
 * @param rad radians
 */
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
	u = sin((lat2r - lat1r) / 2);
	v = sin((lon2r - lon1r) / 2);
	return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

/**
 * Returns the Euclidean distance between two points.
 * @param l1 location in 2D or 3D
 * @param l2 location in 2D or 3D
 * @param l1_index location index
 * @param l2_index location index
 * @param z_flag  0--->2D or 1-->3D 
 * @return The Euclidean distance between the two points
 */
static double calculateDistance(location* l1, location* l2, int l1_index,
		int l2_index, int distance_metric, int z_flag) {

	double z1, z2;
	double x1 = l1->x[l1_index];
	double y1 = l1->y[l1_index];
	double x2 = l2->x[l2_index];
	double y2 = l2->y[l2_index];
	if (l1->z == NULL || l2->z == NULL || z_flag == 1) {
		if (distance_metric == 1)
			return distanceEarth(x1, y1, x2, y2);
		return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
	} else {
		if (distance_metric == 1) {
			printf("Great Circle (GC) distance is only valid for 2d\n");
			exit(0);
		}
		z1 = l1->z[l1_index];
		z2 = l2->z[l2_index];
		return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2));
	}
}

/**
 * Function for spatial range
 * @param x x co-ordinate
 * @param y y co-ordinate
 * @param a parameter for function
 * @param b parameter for function
 * @return The function ae^(sin bx + sin by)
 */
static double lambda(double x, double y, double a, double b) {

	return (a * pow(pow_e, sin(b * x) + sin(b * y)));
}

/**
 * Function for partial sill
 * @param x x co-ordinate
 * @param y y co-ordinate
 * @param d parameter for function
 * @param e parameter for function
 * @param f parameter for function
 * @return The function de^(e(x+y)) + f
 */
static double sigma(double x, double y, double d, double e, double f) {

	return (d * pow(pow_e, e * (x + y)) + f);
}

/**
 * Function for smoothness parameter
 * @param x x co-ordinate
 * @param y y co-ordinate
 * @param g parameter for function
 * @param h parameter for function
 * @param ti parameter for function
 * @return The function ge^(h(x+y)) + ti
 */
static double neu(double x, double y, double g, double h, double ti) {

	return (g * pow(pow_e, h * (x + y)) + ti);
}

/**
 * Utility function that evaluates the matern. Similiar to
 * (https://www.rdocumentation.org/packages/PrevMap/versions/1.5.3/topics/matern.kernel) in R
 * @param range Spatial Range parameter (Also known as rho)
 * @param smoothness Smoothness parameter (Also known as neu)
 * @param distance Distance between the two locations
 * @return Matern function evaluation
 */
static double matern_util(double range, double smoothness, double distance) {
	double con = 0.0;
	con = pow(2, (smoothness - 1)) * tgamma(smoothness);
	con = 1.0 / con;

	if (distance == 0)
		return 1;
	else
		return con * pow(distance / range, smoothness)
			* gsl_sf_bessel_Knu(smoothness, distance / range); // Matern Function
}

/**
 * Returns the Mahalanobis distance between two points to account for anisotropy
 * @param x1 x co-ordinate of first point
 * @param y1 y co-ordinate of first point
 * @param x2 x co-ordinate of second point
 * @param y2 y co-ordinate of second point
 * @param a11 First element of the positive definite matrix that defines the Mahalanobis Distance
 * @param a12 Second element of the positive definite matrix that defines the Mahalanobis Distance
 * @param a21 Third element of the positive definite matrix that defines the Mahalanobis Distance
 * @param a22 Fourth element of the positive definite matrix that defines the Mahalanobis Distance
 * @return The Mahalanobis Distance
 */
static double calculateMahalanobisDistanceSquared(double x1, double y1, double x2, 
		double y2, double a11, double a12, 
		double a21, double a22) {

	double diffx = x1 - x2;
	double diffy = y1 - y2;

	double el1 = a11 * diffx + a21 * diffy;
	double el2 = a12 * diffx + a22 * diffy;

	double ans = el1 * diffx + el2 * diffy;

	return ans;
}


/**
 * Function for TODO
 * @param nu  TODO
 * @param input TODO
 * @return TODO
 */
static double dbessel_input(double nu, double input) {
	return (nu / input * gsl_sf_bessel_Knu(nu, input) - gsl_sf_bessel_Knu(nu + 1, input));
}


/**
 * Function for TODO
 * @param x  TODO
 * @return TODO
 */
static double psi(double x) {
	double result = 0, xx, xx2, xx4;
	assert(x > 0);
	for (; x < 7; ++x)
		result -= 1 / x;
	x -= 1.0 / 2.0;
	xx = 1.0 / x;
	xx2 = xx * xx;
	xx4 = xx2 * xx2;
	result += log(x) + (1. / 24.) * xx2 - (7.0 / 960.0) * xx4 + (31.0 / 8064.0) * xx4 * xx2 -
		(127.0 / 30720.0) * xx4 * xx4;
	return result;
}

/**
 * Function for TODO
 * @param x  TODO
 * @return TODO
 */
static double dpsi(double x) {
	return (1 / x + 0.5 * pow(x, 2));
}

/**
 * Function for TODO
 * @param nu  TODO
 * @param input  TODO
 * @return TODO
 */
static double dbessel_nu(double nu, double input) {
	if (nu == 0) return (0);

	else {
		double dbessel = (gsl_sf_bessel_Knu(nu + 0.000000001, input) - gsl_sf_bessel_Knu(nu, input)) / 0.000000001;
		return (dbessel);
	}
}

/**
 * Function for TODO
 * @param nu  TODO
 * @param input  TODO
 * @return TODO
 */
static double ddbessel_input(double nu, double input) {
	return (-0.5 * (dbessel_input(nu - 1, input) + dbessel_input(nu + 1, input)));
}

/**
 * Function for TODO
 * @param nu  TODO
 * @param input  TODO
 * @return TODO
 */
static double dbessel_input_nu(double nu, double input) {
	if (nu < 1) {
		double nu_new = abs(nu - 1);

		return (-0.5 * (-dbessel_nu(nu_new, input) + dbessel_nu(abs(nu + 1), input)));
	} else {

		return (-0.5 * (dbessel_nu(nu - 1, input) + dbessel_nu(abs(nu + 1), input)));
	}
}

/**
 * Function for TODO
 * @param nu  TODO
 * @param input  TODO
 * @return TODO
 */
static double ddbessel_nu(double nu, double input) {
	if (nu == 0) return (0);
	else {
		double ddbessel = (dbessel_nu(nu + 0.000001, input) - dbessel_nu(nu, input)) / 0.000001;
		return (ddbessel);
	}
}


/*****************************************************************************
 *
 *  core_dcmg_matern_dsigma_square - TODO
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
 *******************************************************************************/

void core_dcmg_matern_dsigma_square(double* A, int m, int n, 
		int m0, int n0, location* l1, 
		location* l2, double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double sigma_square = localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {
				A[i + j * m] = 1 ;
			} else {
				A[i + j * m] = con * pow(expr, localtheta[2]) *
					gsl_sf_bessel_Knu(localtheta[2], expr); // derivative with respect to sigma square
			}

			j0++;
		}
		i0++;
	}
}


/*****************************************************************************
 *
 *  core_dcmg_matern_dnu - TODO
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
 *******************************************************************************/
void core_dcmg_matern_dnu(double* A, int m, int n, 
		int m0, int n0, location* l1, 
		location* l2, double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double nu_expr = 0.0;
	double sigma_square = localtheta[0];
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {
				A[i + j * m] = 0.0;
			} else {

				nu_expr = -2 * log(2.0) * pow(2, -localtheta[2]) * 1 / tgamma(localtheta[2])
					* pow(expr, localtheta[2]) * gsl_sf_bessel_Knu(localtheta[2], expr) +
					pow(2, 1 - localtheta[2])
					* (-1 / tgamma(localtheta[2]) * psi(localtheta[2]) * pow(expr, localtheta[2])
							* gsl_sf_bessel_Knu(localtheta[2], expr) + 1 / tgamma(localtheta[2])
							* (pow(expr, localtheta[2]) * log(expr)
								* gsl_sf_bessel_Knu(localtheta[2], expr) +
								pow(expr, localtheta[2])
								* dbessel_nu(localtheta[2], expr)));

				A[i + j * m] = sigma_square * nu_expr; //derivative with respect to nu

			}
			j0++;

		}
		i0++;
	}
}

/*****************************************************************************
 *
 *  core_dcmg_matern_dbeta - TODO
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
 *******************************************************************************/
void core_dcmg_matern_dbeta(double* A, int m, int n, 
		int m0, int n0, location* l1, 
		location* l2, double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double beta_expr = 0.0;
	double con = 0.0;
	double sigma_square = localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;

	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {
				A[i + j * m] = 0.0; 
			} else {
				beta_expr = -localtheta[2] / localtheta[1] * pow(expr, localtheta[2])
					* gsl_sf_bessel_Knu(localtheta[2], expr) - pow(expr, localtheta[2])
					* dbessel_input(localtheta[2], expr) * expr /
					localtheta[1];

				A[i + j * m] = sigma_square * con * beta_expr;//derivative withrespect to beta

			}
			j0++;
		}


		i0++;
	}
}


/*****************************************************************************
 *
 *  core_dcmg_matern_ddsigma_square - TODO
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
 *******************************************************************************/
void core_dcmg_matern_ddsigma_square(double* A, int m, int n) {
	int i, j;
	double x0, y0, z0;

	for (i = 0; i < m; i++) {

		for (j = 0; j < n; j++) {

			A[i + j * m] = 0.0;

		}
	}
}

/*****************************************************************************
 *
 *  core_dcmg_matern_ddsigma_square_beta - TODO
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
 *******************************************************************************/
void core_dcmg_matern_ddsigma_square_beta(double* A, int m, int n, 
		int m0, int n0, location* l1,
		location* l2,	double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double beta_expr = 0.0;
	double sigma_square = localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {

				A[i + j * m] = 0.0;

			} else {
				beta_expr = -localtheta[2] / localtheta[1] * pow(expr, localtheta[2])
					* gsl_sf_bessel_Knu(localtheta[2], expr) - pow(expr, localtheta[2])
					* dbessel_input(localtheta[2], expr) * expr /
					localtheta[1];
				A[i + j * m] = con * beta_expr; // derivative with respect to sigma square and beta
			}
			j0++;
		}
		i0++;
	}
}

/*****************************************************************************
 *
 *  core_dcmg_matern_ddsigma_square_nu - TODO
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
 *******************************************************************************/
void core_dcmg_matern_ddsigma_square_nu(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2,	double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double nu_expr = 0.0;
	double sigma_square = localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {
				A[i + j * m] = 0.0;
			} else {

				nu_expr = (1 - localtheta[2]) * 1 / pow(2, localtheta[2]) * 1 / tgamma(localtheta[2])
					* pow(expr, localtheta[2]) * gsl_sf_bessel_Knu(localtheta[2], expr) +
					pow(2, 1 - localtheta[2])
					* (-1 / tgamma(localtheta[2]) * psi(localtheta[2]) * pow(expr, localtheta[2])
							* gsl_sf_bessel_Knu(localtheta[2], expr) + 1 / tgamma(localtheta[2])
							* (pow(expr, localtheta[2]) * log(expr)
								* gsl_sf_bessel_Knu(localtheta[2], expr) +
								pow(expr, localtheta[2])
								* dbessel_nu(localtheta[2], expr)));
				A[i + j * m] = nu_expr; // derivative with respect to sigma square
			}
			j0++;
		}
		i0++;
	}
}


/*****************************************************************************
 *
 *  core_dcmg_matern_ddbeta_beta - TODO
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
 *******************************************************************************/
void core_dcmg_matern_ddbeta_beta(double* A, int m, int n, 
		int m0, int n0, location* l1, 
		location* l2, double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double beta_expr = 0.0;
	double beta_expr_prime = 0.0;
	double sigma_square = localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {

				A[i + j * m] = 0.0;
			} else {

				beta_expr = -localtheta[2] / localtheta[1] * pow(expr, localtheta[2])
					* gsl_sf_bessel_Knu(localtheta[2], expr) - pow(expr, localtheta[2])
					* dbessel_input(localtheta[2], expr) * expr /
					localtheta[1];

				beta_expr_prime = -localtheta[2] / localtheta[1] * pow(expr, localtheta[2])
					* dbessel_input(localtheta[2], expr) - pow(expr, localtheta[2])
					* ddbessel_input(localtheta[2], expr) * expr /
					localtheta[1];

				A[i + j * m] = (localtheta[2] / pow(localtheta[1], 2) * pow(expr, localtheta[2]) *
						gsl_sf_bessel_Knu(localtheta[2], expr)
						- localtheta[2] / localtheta[1] * beta_expr +
						2 * expr / pow(localtheta[1], 2) * pow(expr, localtheta[2])
						* dbessel_input(localtheta[2], expr) - expr / localtheta[1] * beta_expr_prime) *
					sigma_square * con; // derivative with respect to beta beta
			}
			j0++;
		}
		i0++;
	}
}


/*****************************************************************************
 *
 *  core_dcmg_matern_ddbeta_beta - TODO
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
 *******************************************************************************/
void core_dcmg_matern_ddbeta_nu(double* A, int m, int n, 
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double nu_expr = 0.0;
	double nu_expr_prime = 0.0;
	double sigma_square = localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {
				A[i + j * m] = 0.0;
			} else {
				nu_expr = (1 - localtheta[2]) * 1 / pow(2, localtheta[2]) * 1 / tgamma(localtheta[2])
					* pow(expr, localtheta[2]) * gsl_sf_bessel_Knu(localtheta[2], expr) +
					pow(2, 1 - localtheta[2])
					* (-1 / tgamma(localtheta[2]) * psi(localtheta[2]) * pow(expr, localtheta[2])
							* gsl_sf_bessel_Knu(localtheta[2], expr) + 1 / tgamma(localtheta[2])
							* (pow(expr, localtheta[2]) * log(expr)
								* gsl_sf_bessel_Knu(localtheta[2], expr) +
								pow(expr, localtheta[2])
								* dbessel_nu(localtheta[2], expr)));

				nu_expr_prime = (1 - localtheta[2]) * 1 / pow(2, localtheta[2]) * 1 / tgamma(localtheta[2])
					* pow(expr, localtheta[2]) * dbessel_input(localtheta[2], expr) +
					pow(2, 1 - localtheta[2])
					* (-1 / tgamma(localtheta[2]) * psi(localtheta[2]) * pow(expr, localtheta[2])
							* dbessel_input(localtheta[2], expr) + 1 / tgamma(localtheta[2])
							* (pow(expr, localtheta[2]) * log(expr)
								* dbessel_input(localtheta[2], expr) +
								pow(expr, localtheta[2])
								* dbessel_input_nu(localtheta[2], expr)));

				A[i + j * m] = (-1 / localtheta[1] * (con * pow(expr, localtheta[2])
							* gsl_sf_bessel_Knu(localtheta[2], expr))
						- localtheta[2] / localtheta[1] * nu_expr - expr / localtheta[1] * nu_expr_prime) *
					sigma_square;
			}
			j0++;
		}
		i0++;
	}
}

/*****************************************************************************
 *
 *  core_dcmg_matern_ddnu_nu - TODO
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
 *******************************************************************************/
void core_dcmg_matern_ddnu_nu(double* A, int m, int n, 
		int m0, int n0, location* l1, 
		location* l2, double* localtheta, int distance_metric) {
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double nu_expr = 0.0;
	double nu_expr_dprime = 0.0;
	double sigma_square = localtheta[0];//*  localtheta[0];
	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0) {

				A[i + j * m] = 0.0 /*+ 1e-4*/;
			} else {
				nu_expr = (1 - localtheta[2]) * 1 / pow(2, localtheta[2]) * 1 / tgamma(localtheta[2])
					* pow(expr, localtheta[2]) * gsl_sf_bessel_Knu(localtheta[2], expr) +
					pow(2, 1 - localtheta[2])
					* (-1 / tgamma(localtheta[2]) * psi(localtheta[2]) * pow(expr, localtheta[2])
							* gsl_sf_bessel_Knu(localtheta[2], expr) + 1 / tgamma(localtheta[2])
							* (pow(expr, localtheta[2]) * log(expr)
								* gsl_sf_bessel_Knu(localtheta[2], expr) +
								pow(expr, localtheta[2])
								* dbessel_nu(localtheta[2], expr)));


				nu_expr_dprime = (1 - localtheta[2]) * 1 / pow(2, localtheta[2]) * 1 / tgamma(localtheta[2])
					* pow(expr, localtheta[2]) * dbessel_nu(localtheta[2], expr) +
					pow(2, 1 - localtheta[2])
					* (-1 / tgamma(localtheta[2]) * psi(localtheta[2]) * pow(expr, localtheta[2])
							* dbessel_nu(localtheta[2], expr) + 1 / tgamma(localtheta[2])
							* (pow(expr, localtheta[2]) * log(expr)
								* dbessel_nu(localtheta[2], expr) +
								pow(expr, localtheta[2])
								* ddbessel_nu(localtheta[2], expr)));

				A[i + j * m] = (-0.5 * con * pow(expr, localtheta[2]) * gsl_sf_bessel_Knu(localtheta[2], expr) +
						(1 - localtheta[2]) / 2 * nu_expr
						- dpsi(localtheta[2]) * con * pow(expr, localtheta[2]) *
						gsl_sf_bessel_Knu(localtheta[2], expr) - psi(localtheta[2])
						* nu_expr + log(expr) * nu_expr +
						nu_expr_dprime) * sigma_square;
				// derivative with respect to nu nu
			}
			j0++;
		}
		i0++;

	}

}


/*****************************************************************************
 *
 *  core_dcmg - TODO
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
 *******************************************************************************/
void core_dcmg(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double sigma_square = localtheta[0];// * localtheta[0];

	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	con = sigma_square * con;

	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0)
				A[i + j * m] = sigma_square /*+ 1e-4*/;
			else
				A[i + j * m] = con * pow(expr, localtheta[2])
					* gsl_sf_bessel_Knu(localtheta[2], expr); // Matern Function

			j0++;
		}
		i0++;
	}

}

/*****************************************************************************
 *
 *  core_ng_dcmg - TODO
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
 *******************************************************************************/
void core_ng_dcmg(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	//localtheta[0] <- \phi
	//localtheta[1] <- \nu
	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double sigma_square = 1;

	con = pow(2, (localtheta[1] - 1)) * tgamma(localtheta[1]);
	con = 1.0 / con;
	con = sigma_square * con;

	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = 4 * sqrt(2 * localtheta[1]) *
				(calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[0]);
			if (expr == 0)
				A[i + j * m] = sigma_square /*+ 1e-4*/;
			else
				A[i + j * m] = con * pow(expr, localtheta[1])
					* gsl_sf_bessel_Knu(localtheta[1], expr); // Matern Function
			j0++;
		}
		i0++;
	}
}

/*****************************************************************************
 *
 *  core_ng_exp_dcmg - TODO
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
 *******************************************************************************/
void core_ng_exp_dcmg(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	//localtheta[0] <- \phi
	//localtheta[1] <- \nu

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;

	double sigma_square = 1;

	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[0];

			if (expr == 0)
				A[i + j * m] = sigma_square /*+ 1e-4*/;
			else
				A[i + j * m] = exp(-expr);

			j0++;
		}
		i0++;
	}
}

/***************************************************************************//**
 *
 *  core_scmg - TODO
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
void core_scmg(float *A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	float x0, y0, z0;
	float expr = 0.0;
	float con = 0.0;
	float sigma_square = localtheta[0];

	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	con = sigma_square * con;


	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0)
				A[i + j * m] = sigma_square /*+ 1e-4*/;
			else
				A[i + j * m] = con * pow(expr, localtheta[2])
					* gsl_sf_bessel_Knu(localtheta[2], expr); // Matern Function
			j0++;
		}
		i0++;
	}
}

/*****************************************************************************
 *
 *  core_sdcmg - TODO
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
 *******************************************************************************/
void core_sdcmg(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0;
	double expr = 0.0;
	double con = 0.0;
	double sigma_square = localtheta[0];

	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	con = sigma_square * con;


	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[1];
			if (expr == 0)
				A[i + j * m] = (float) (sigma_square /*+ 1e-4*/);
			else
				A[i + j * m] = (float) (con * pow(expr, localtheta[2])
						* gsl_sf_bessel_Knu(localtheta[2], expr)); // Matern Function
			j0++;
		}
		i0++;
	}
}

/*******************************************************************************
 *
 *  core_dcmg_spacetime_matern - TODO.
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
 *          localtheta[0] is the variance of Z.
 *          localtheta[1] is the range parameter in space of Z.
 *          localtheta[2] is the smoothness parameter in space of Z.
 *          localtheta[3] is the range parameter in time of Z.
 *          localtheta[4] is the smoothness parameter in time of Z.
 *          localtheta[5] is the space-time nonseparability parameter of Z.
 *          localtheta[6] is an auxilliary parameter in time of Z.
 *
 *          localtheta[6] is an auxilliary parameter in time of Z.
 *
 *          Reference: Equation (16) of Tilmann Gneiting (2002) Nonseparable, 
 * Stationary Covariance Functions for Space–Time Data, 
 * Journal of the American Statistical Association, 97:458, 590-600, DOI: 10.1198/016214502760047113
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
void core_dcmg_spacetime_matern(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0, z1;
	double expr = 0.0, expr1 = 0.0, expr2 = 0.0, expr3 = 0.0, expr4 = 0.0;
	double con = 0.0;
	double sigma_square = localtheta[0];

	con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
	con = 1.0 / con;
	con = sigma_square * con;

	for (i = 0; i < m; i++) {
		j0 = n0;
		z0 = l1->z[i0];
		for (j = 0; j < n; j++) {
			z1 = l2->z[j0];

			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 1) / localtheta[1];
			expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * localtheta[4]) / localtheta[3] + 1.0, localtheta[5] / 2.0);
			expr3 = expr / expr2;
			expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * localtheta[4]) / localtheta[3] + 1.0,
					localtheta[5] + localtheta[6]);

			if (expr == 0)
				A[i + j * m] = sigma_square / expr4 /*+ 1e-4*/;
			else
				A[i + j * m] = con * pow(expr3, localtheta[2])
					* gsl_sf_bessel_Knu(localtheta[2], expr3) / expr4; // Matern Function
			j0++;
		}
		i0++;
	}
}

/*****************************************************************************
 *
 *  core_dcmg_spacetime_bivariate_parsimonious - TODO
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
 *******************************************************************************/
void core_dcmg_spacetime_bivariate_parsimonious(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0, z0, z1;
	double expr = 0.0, expr1 = 0.0, expr2 = 0.0, expr3 = 0.0, expr4 = 0.0;
	double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;

	con1 = pow(2, (localtheta[3] - 1)) * tgamma(localtheta[3]);
	con1 = 1.0 / con1;
	con1 = localtheta[0] * con1;

	con2 = pow(2, (localtheta[4] - 1)) * tgamma(localtheta[4]);
	con2 = 1.0 / con2;
	con2 = localtheta[1] * con2;

	//The average
	nu12 = 0.5 * (localtheta[3] + localtheta[4]);

	rho = localtheta[5] * sqrt((tgamma(localtheta[3] + 1) * tgamma(localtheta[4] + 1)) /
			(tgamma(localtheta[3]) * tgamma(localtheta[4]))) *
		tgamma(nu12) / tgamma(nu12 + 1);

	con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
	con12 = 1.0 / con12;
	con12 = rho * sqrt(localtheta[0] * localtheta[1]) * con12;

	i0 /= 2;
	for (i = 0; i < m - 1; i += 2) {
		j0 = n0 / 2;
		z0 = l1->z[i0];

		for (j = 0; j < n - 1; j += 2) {
			z1 = l2->z[j0];

			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 1) / (localtheta[2] * 1000);
			expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * localtheta[7]) / localtheta[6] + 1, localtheta[8] / 2);
			expr3 = expr / expr2;
			expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * localtheta[7]) / localtheta[6] + 1,
					localtheta[8] + localtheta[9]);

			if (expr == 0) {
				A[i + j * m] = localtheta[0] / expr4;
				A[(i + 1) + j * m] = A[i + (j + 1) * m] = rho
					* sqrt(localtheta[0] * localtheta[1]) / expr4;
				A[(i + 1) + (j + 1) * m] = localtheta[1] / expr4;
			} else {
				A[i + j * m] = con1 * pow(expr3, localtheta[3])
					* gsl_sf_bessel_Knu(localtheta[3], expr3) / expr4;
				A[(i + 1) + j * m] = A[i + (j + 1) * m] = con12 * pow(expr3, nu12)
					* gsl_sf_bessel_Knu(nu12, expr3) / expr4;
				A[(i + 1) + (j + 1) * m] = con2 * pow(expr3, localtheta[4])
					* gsl_sf_bessel_Knu(localtheta[4], expr3) / expr4;
			}
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
 *          localtheta[0] is the scale of Z1.
 *          localtheta[1] is the scale of Z2.
 *          localtheta[2] is the first latent cross scale parameter; see Delta_B of Remark 1 (c) of Apanasovich et al. (2012).
 *          localtheta[3] is the second latent cross scale parameter; see B_ij of Remark 1 (c) of Apanasovich et al. (2012).
 *          localtheta[4] is the smoothness of Z1.
 *          localtheta[5] is the smoothness of Z2.
 *          localtheta[6] is the first latent cross smoothness parameter; see Delta_A of Theorem 1 (i) of Apanasovich et al. (2012).
 *          localtheta[7] is the second latent cross smoothness parameter; see A_ij of Theorem 1 (i) of Apanasovich et al. (2012).
 *          localtheta[8] is the first latent cross variance parameter; see W_ii of Equation (8) of Apanasovich et al. (2012).
 *          localtheta[9] is the second latent cross variance parameter; see W_jj of Equation (8) of Apanasovich et al. (2012).
 *          localtheta[10] is the third latent cross variance parameter; see W_{V,ij} of Equation (8) of Apanasovich et al. (2012).
 *
 *
 *          Reference: Tatiyana V. Apanasovich , Marc G. Genton & Ying Sun (2012) A Valid Matérn Class of Cross-Covariance Functions 
 *          for Multivariate Random Fields With Any Number of Components, Journal of the American Statistical Association,
 *          107:497, 180-193, DOI: 10.1080/01621459.2011.643197
 *
 * @param[in] distance_metric
 *          Distance metric "euclidean Distance (ED) ->0" or "Great Circle Distance (GCD) ->1"
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/

void core_dcmg_bivariate_flexible(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0;
	double expr1 = 0.0, expr2 = 0.0, expr12 = 0.0;
	double con1 = 0.0, con2 = 0.0, con12 = 0.0, scale12 = 0.0, rho = 0.0, nu12 = 0.0, sigma_square11 = 0.0, sigma_square22 = 0.0;
	double scale1 = localtheta[0], scale2 = localtheta[1], nu1 = localtheta[4], nu2 = localtheta[5];

	scale12 = pow(0.5 * (pow(scale1, -2) + pow(scale2, -2)) + localtheta[2] * (1 - localtheta[3]),
			-0.5); //Remark 1 (c) of Apanasovich et al. (2012)

	nu12 = 0.5 * (nu1 + nu2) + localtheta[6] * (1 - localtheta[7]); //Theorem 1 (i) of Apanasovich et al. (2012).

	rho = localtheta[8] * localtheta[9] * localtheta[10] *
		pow(scale12, 2 * localtheta[6] + (nu1 + nu2))
		* tgamma(0.5 * (nu1 + nu2) + 1) * tgamma(nu12) /
		tgamma(nu12 + 1); //Equation (8) of Apanasovich et al. (2012).

	sigma_square11 = localtheta[8] * localtheta[8] *
		pow(scale1, 2 * localtheta[6] + nu1 + nu1) *
		tgamma(nu1); //Equation (8) of Apanasovich et al. (2012).

	sigma_square22 = localtheta[9] * localtheta[9] *
		pow(scale2, 2 * localtheta[6] + nu2 + nu2) *
		tgamma(nu2); //Equation (8) of Apanasovich et al. (2012).

	con1 = pow(2, (nu1 - 1)) * tgamma(nu1);
	con1 = 1.0 / con1;
	con1 = sigma_square11 * con1;

	con2 = pow(2, (nu2 - 1)) * tgamma(nu2);
	con2 = 1.0 / con2;
	con2 = sigma_square22 * con2;

	con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
	con12 = 1.0 / con12;
	con12 = rho * con12;

	i0 /= 2;
	for (i = 0; i < m; i += 2) {
		j0 = n0 / 2;
		for (j = 0; j < n; j += 2) {
			expr1 = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / scale1;
			expr2 = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / scale2;
			expr12 = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / scale12;

			if (expr1 == 0) {
				A[i + j * m] = localtheta[0];
				A[(i + 1) + j * m] = A[i + (j + 1) * m] = rho;
				A[(i + 1) + (j + 1) * m] = localtheta[1];
			} else {
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
void core_dcmg_bivariate_parsimonious(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0;
	double expr = 0.0;
	double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;

	con1 = pow(2, (localtheta[3] - 1)) * tgamma(localtheta[3]);
	con1 = 1.0 / con1;
	con1 = localtheta[0] * con1;

	con2 = pow(2, (localtheta[4] - 1)) * tgamma(localtheta[4]);
	con2 = 1.0 / con2;
	con2 = localtheta[1] * con2;

	//The average
	nu12 = 0.5 * (localtheta[3] + localtheta[4]);

	rho = localtheta[5] * sqrt((tgamma(localtheta[3] + 1) * tgamma(localtheta[4] + 1)) /
			(tgamma(localtheta[3]) * tgamma(localtheta[4]))) *
		tgamma(nu12) / tgamma(nu12 + 1);


	con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
	con12 = 1.0 / con12;
	con12 = rho * sqrt(localtheta[0] * localtheta[1]) * con12;

	i0 /= 2;
	for (i = 0; i < m - 1; i += 2) {
		j0 = n0 / 2;
		for (j = 0; j < n - 1; j += 2) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[2];

			if (expr == 0) {
				A[i + j * m] = localtheta[0];
				A[(i + 1) + j * m] = A[i + (j + 1) * m] = rho
					* sqrt(localtheta[0] * localtheta[1]);
				A[(i + 1) + (j + 1) * m] = localtheta[1];
			} else {
				A[i + j * m] = con1 * pow(expr, localtheta[3])
					* gsl_sf_bessel_Knu(localtheta[3], expr);
				A[(i + 1) + j * m] = A[i + (j + 1) * m] = con12 * pow(expr, nu12)
					* gsl_sf_bessel_Knu(nu12, expr);
				A[(i + 1) + (j + 1) * m] = con2 * pow(expr, localtheta[4])
					* gsl_sf_bessel_Knu(localtheta[4], expr);
			}
			j0++;
		}
		i0++;
	}
}


/***************************************************************************//**
 *
 *  core_dcmg_bivariate_parsimonious - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Parsimonious trivariate Matern Kernel).
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
void core_dcmg_trivariate_parsimonious(double* A, int m, int n,
		int m0, int n0, location* l1,
		location* l2, double* localtheta, int distance_metric) {

	int i, j;
	int i0 = m0;
	int j0 = n0;
	double x0, y0;
	double expr = 0.0;
	double con1 = 0.0, con2 = 0.0, con3 = 0.0, con12 = 0.0,
		   con13 = 0.0, con23 = 0.0, rho12 = 0.0, rho13 = 0.0, rho23 = 0.0, nu12 = 0.0, nu13 = 0.0, nu23 = 0.0;

	con1 = pow(2, (localtheta[4] - 1)) * tgamma(localtheta[4]);
	con1 = 1.0 / con1;
	con1 = localtheta[0] * con1;

	con2 = pow(2, (localtheta[5] - 1)) * tgamma(localtheta[5]);
	con2 = 1.0 / con2;
	con2 = localtheta[1] * con2;

	con3 = pow(2, (localtheta[6] - 1)) * tgamma(localtheta[6]);
	con3 = 1.0 / con3;
	con3 = localtheta[2] * con3;

	//The average
	nu12 = 0.5 * (localtheta[4] + localtheta[5]);
	nu13 = 0.5 * (localtheta[4] + localtheta[6]);
	nu23 = 0.5 * (localtheta[5] + localtheta[6]);

	rho12 = localtheta[7] * sqrt((tgamma(localtheta[4] + 1) * tgamma(localtheta[5] + 1)) /
			(tgamma(localtheta[4]) * tgamma(localtheta[5]))) *
		tgamma(nu12) / tgamma(nu12 + 1);

	rho13 = localtheta[8] * sqrt((tgamma(localtheta[4] + 1) * tgamma(localtheta[6] + 1)) /
			(tgamma(localtheta[4]) * tgamma(localtheta[6]))) *
		tgamma(nu13) / tgamma(nu13 + 1);

	rho23 = localtheta[9] * sqrt((tgamma(localtheta[5] + 1) * tgamma(localtheta[6] + 1)) /
			(tgamma(localtheta[5]) * tgamma(localtheta[6]))) *
		tgamma(nu23) / tgamma(nu23 + 1);

	con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
	con12 = 1.0 / con12;
	con12 = rho12 * sqrt(localtheta[0] * localtheta[1]) * con12;

	con13 = pow(2, (nu13 - 1)) * tgamma(nu13);
	con13 = 1.0 / con13;
	con13 = rho13 * sqrt(localtheta[0] * localtheta[2]) * con13;

	con23 = pow(2, (nu23 - 1)) * tgamma(nu23);
	con23 = 1.0 / con23;
	con23 = rho23 * sqrt(localtheta[1] * localtheta[2]) * con23;

	i0 /= 3;
	for (i = 0; i < m - 1; i += 3) {
		j0 = n0 / 3;
		for (j = 0; j < n - 1; j += 3) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0) / localtheta[3];

			if (expr == 0) {
				A[i + j * m] = localtheta[0];

				A[(i + 1) + j * m] = A[i + (j + 1) * m] = rho12 * sqrt(localtheta[0] * localtheta[1]);
				A[(i + 2) + j * m] = A[i + (j + 2) * m] = rho13 * sqrt(localtheta[0] * localtheta[2]);

				A[(i + 1) + (j + 1) * m] = localtheta[1];

				A[(i + 1) + (j + 2) * m] = A[(i + 2) + (j + 1) * m] = rho23 * sqrt(localtheta[1] * localtheta[2]);

				A[(i + 2) + (j + 2) * m] = localtheta[2];
			} else {
				A[i + j * m] = con1 * pow(expr, localtheta[4]) * gsl_sf_bessel_Knu(localtheta[4], expr);

				A[(i + 1) + j * m] = A[i + (j + 1) * m] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12, expr);
				A[(i + 2) + j * m] = A[i + (j + 2) * m] = con13 * pow(expr, nu13) * gsl_sf_bessel_Knu(nu13, expr);

				A[(i + 1) + (j + 1) * m] = con2 * pow(expr, localtheta[5]) * gsl_sf_bessel_Knu(localtheta[5], expr);

				A[(i + 1) + (j + 2) * m] = A[(i + 2) + (j + 1) * m] =
					con23 * pow(expr, nu23) * gsl_sf_bessel_Knu(nu23, expr);

				A[(i + 2) + (j + 2) * m] = con3 * pow(expr, localtheta[6]) * gsl_sf_bessel_Knu(localtheta[6], expr);
			}
			j0++;
		}
		i0++;
	}
}

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
void core_dcmg_pow_exp(double* A, int m, int n,
		int m0, int n0,
		location* l1, location* l2,
		double* localtheta, int distance_metric) {

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
	double expr = 0.0;
	double expr1 = 0.0;
	double sigma_square = localtheta[0];// * localtheta[0];

	for (i = 0; i < m; i++) {
		j0 = n0;
		for (j = 0; j < n; j++) {
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0);
			expr1 = pow(expr, localtheta[2]);


			if (expr == 0)
				A[i + j * m] = sigma_square /*+ 1e-4*/;
			else
				A[i + j * m] = sigma_square * exp(-(expr1 / localtheta[1])); // power-exp kernel
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
void core_scmg_pow_exp(float *A, int m, int n,
		int m0, int n0,
		location* l1, location* l2,
		double* localtheta, int distance_metric) {

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
			expr = calculateDistance(l1, l2,
					i0, j0, distance_metric, 0);
			expr1 = pow(expr, localtheta[2]);
			if (expr == 0)
				A[i + j * m] = sigma_square /*+ 1e-4*/;
			else
				A[i + j * m] = sigma_square * exp(-(expr1 / localtheta[1])); // power-exp kernel
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

void core_sdcmg_pow_exp(double* A, int m, int n,
		int m0, int n0,
		location* l1, location* l2,
		double* localtheta, int distance_metric) {

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
			expr = calculateDistance(l1, l2, i0, j0, distance_metric, 0);
			expr1 = pow(expr, localtheta[2]);
			if (expr == 0)
				A[i + j * m] = (float) (sigma_square /*+ 1e-4*/);
			else
				A[i + j * m] = (float) (sigma_square * exp(-(expr1 / localtheta[1]))); // power-exp kernel
			j0++;
		}
		i0++;
	}
}

/***************************************************************************//**
 *
 *  core_dcmg_non_stat - Generate covariance matrix A in dense format between two sets of locations (l1, l2) (Non-Stationary Matern Kernel).
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
void core_dcmg_non_stat(double* A, int m, int n, int m0, int n0, location* l1, location* l2, double* localtheta,
		int distance_metric) {

	double l1x, l1y, l2x, l2y;

	double expr = 0.0;
	double con, sigma_square, beta, nu;
	double a, b, c, d, e, f, g, h, ti, tj;

	a = localtheta[0];
	b = localtheta[1];
	d = localtheta[2];
	e = localtheta[3];
	f = localtheta[4];
	g = localtheta[5];
	h = localtheta[6];
	ti = localtheta[7];
	tj = 0;

	double nu_arr_1[m];
	double sigma_arr_1[m];
	double lambda_arr_1[m];
	double nu_arr_2[n];
	double sigma_arr_2[n];
	double lambda_arr_2[n];
	for (int i = 0; i < m; i++) {
		l1x = l1->x[i + m0];
		l1y = l1->y[i + m0];
		nu_arr_1[i] = neu(l1x, l1y, g, h, ti);
		sigma_arr_1[i] = sigma(l1x, l1y, d, e, f);
		lambda_arr_1[i] = lambda(l1x, l1y, a, b);
	}

	for (int j = 0; j < n; j++) {
		l2x = l2->x[j + n0];
		l2y = l2->y[j + n0];

		nu_arr_2[j] = neu(l2x, l2y, g, h, ti);
		sigma_arr_2[j] = sigma(l2x, l2y, d, e, f);
		lambda_arr_2[j] = lambda(l2x, l2y, a, b);
	}

	for (int i = 0; i < m; i++) {
		l1x = l1->x[i + m0];
		l1y = l1->y[i + m0];

		for (int j = 0; j < n; j++) {
			l2x = l2->x[j + n0];
			l2y = l2->y[j + n0];

			double term1 = (sigma_arr_1[i]) * (sigma_arr_2[j]) * sqrt(lambda_arr_1[i]) * sqrt(lambda_arr_2[j]);
			double term2 = 2 / ((lambda_arr_1[i]) + (lambda_arr_2[j]));
			double neuij = ((nu_arr_1[i]) + (nu_arr_2[j])) / 2;
			double Qij = calculateMahalanobisDistanceSquared(l1x, l1y, l2x, l2y, term2, 0, 0, term2);
			double prod1 = 2 * sqrt(neuij * Qij);
			double term3 = matern_util(1, neuij, prod1);
			A[i + j * m] = term1 * term2 * term3;
		}
	}
}
