/**
 *
 * @file testing_dposv.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Bilel Hadri, Hatem Ltaief
 * @date 2010-11-15
 * @generated d Fri Jan 22 15:11:13 2016
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <dlib/optimization.h>
#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "flops.h"

//Sameh Testing MOde
#define MIN_RAND       -0.4
#define MAX_RAND        0.4
#define R        2
#define L       3
#define THETA 0.1//weak=0.03, medium=0.1,  strong=0.3
#define USAGE(args, details)                       \
  printf(" Proper Usage is : ./main_genton "args" \n" \
         details);

#define PI (3.141592653589793)

#define START_TIMING(_t)                          \
  _t =- cWtime();

#define STOP_TIMING(_t)                           \
  _t += cWtime();

double cWtime(void) {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + 1e-6 * tp.tv_usec;
}

enum blas_order_type {
	blas_rowmajor = 101, blas_colmajor = 102
};

enum blas_cmach_type {
	blas_base = 151,
	blas_t = 152,
	blas_rnd = 153,
	blas_ieee = 154,
	blas_emin = 155,
	blas_emax = 156,
	blas_eps = 157,
	blas_prec = 158,
	blas_underflow = 159,
	blas_overflow = 160,
	blas_sfmin = 161
};

enum blas_norm_type {
	blas_one_norm = 171,
	blas_real_one_norm = 172,
	blas_two_norm = 173,
	blas_frobenius_norm = 174,
	blas_inf_norm = 175,
	blas_real_inf_norm = 176,
	blas_max_norm = 177,
	blas_real_max_norm = 178
};

static void BLAS_error(char *rname, int err, int val, int x) {
	fprintf( stderr, "%s %d %d %d\n", rname, err, val, x);
	abort();
}

static
void BLAS_dge_norm(enum blas_order_type order, enum blas_norm_type norm, int m,
		int n, const double *a, int lda, double *res) {
	int i, j;
	float anorm, v;
	char rname[] = "BLAS_dge_norm";

	if (order != blas_colmajor)
		BLAS_error(rname, -1, order, 0);

	if (norm == blas_frobenius_norm) {
		anorm = 0.0f;
		for (j = n; j; --j) {
			for (i = m; i; --i) {
				v = a[0];
				anorm += v * v;
				a++;
			}
			a += lda - m;
		}
		anorm = sqrt(anorm);
	} else if (norm == blas_inf_norm) {
		anorm = 0.0f;
		for (i = 0; i < m; ++i) {
			v = 0.0f;
			for (j = 0; j < n; ++j) {
				v += fabs(a[i + j * lda]);
			}
			if (v > anorm)
				anorm = v;
		}
	} else {
		BLAS_error(rname, -2, norm, 0);
		return;
	}

	if (res)
		*res = anorm;
}

static
double BLAS_dpow_di(double x, int n) {
	double rv = 1.0;

	if (n < 0) {
		n = -n;
		x = 1.0 / x;
	}

	for (; n; n >>= 1, x *= x) {
		if (n & 1)
			rv *= x;
	}

	return rv;
}

static
double BLAS_dfpinfo(enum blas_cmach_type cmach) {
	double eps = 1.0, r = 1.0, o = 1.0, b = 2.0;
	int t = 53, l = 1024, m = -1021;
	char rname[] = "BLAS_dfpinfo";

	if ((sizeof eps) == sizeof(float)) {
		t = 24;
		l = 128;
		m = -125;
	} else {
		t = 53;
		l = 1024;
		m = -1021;
	}

	/* for (i = 0; i < t; ++i) eps *= half; */
	eps = BLAS_dpow_di(b, -t);
	/* for (i = 0; i >= m; --i) r *= half; */
	r = BLAS_dpow_di(b, m - 1);

	o -= eps;
	/* for (i = 0; i < l; ++i) o *= b; */
	o = (o * BLAS_dpow_di(b, l - 1)) * b;

	switch (cmach) {
	case blas_eps:
		return eps;
	case blas_sfmin:
		return r;
	default:
		BLAS_error(rname, -1, cmach, 0);
		break;
	}
	return 0.0;
}

static int check_factorization(int, double*, double*, int, int, double);
//static int check_solution(int, int, double*, int, double*, double*, int, double);

static int GenerateDistanceMatrix(double* dist_matrix, int n);
static int GenerateCovMatrix(double* C, double * dist_matrix, int N,
		double theta);
static double CalculateDistance(double x1, double y1, double x2, double y2);
static double uniform_distribution(double rangeLow, double rangeHigh);
//static void writetofile(char* fname, int m, int n, double* a, int lda);
//static void print_matrix(char* desc, int m, int n, double* a, int lda);

int PLASMA_ML_GenerateZVec_tile();
double PLASMA_MLE(double theta);
void PLASMA_ML_finalize();


double * C ;
double * dist_matrix;
double * Z;
double * Zcpy;
double * Ccpy;
int N,NRHS,LDC,LDZ,check,verbose;
int iter_count=0;


int main(int argc, char **argv) {

	/* Check for number of arguments*/
	if (argc != 6) {
		USAGE("NCORES N TS CHECK VERBOSE",
				"   - NCORES   : number of cores\n" "   - N        : the size of the matrix\n" "   - TS       : tile size\n" //sameh
				"   - CHECK    : check the factorization and the solution\n" "   - VERBOSE  : verbose\n");
		return -1;
	}

	int ncores = atoi(argv[1]);
	N = atoi(argv[2]);
	int ts = atoi(argv[3]);
        check = atoi(argv[4]);
 	verbose = atoi(argv[5]);
	/*-------------------------------------------------------------
	*		Initialization of PLASMA
	*/
	PLASMA_Init(ncores);
	PLASMA_Disable(PLASMA_AUTOTUNING);
	PLASMA_Set(PLASMA_TILE_SIZE, ts);
	PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING);
	PLASMA_Enable(PLASMA_WARNINGS);
	PLASMA_Enable(PLASMA_ERRORS);
	/*-------------------------------------------------------------*/

	double time_opt=0.0;
	double starting_theta=0.5;
	double max_theta_hat=0.0;
	NRHS = 1;
	LDC = N;
	LDZ = N;

	//Memory Allocation
	C = (double *) malloc(LDC * N * sizeof(double));
	dist_matrix = (double *) malloc(LDC * N * sizeof(double));
	Z = (double *) malloc(LDZ * NRHS * sizeof(double));
	Zcpy = (double *) malloc(LDZ * NRHS * sizeof(double));


	if (check == 1)
		Ccpy = (double *) malloc(LDC * N * sizeof(double));

	/* Check if unable to allocate memory */
	if ((!C) || (!Z) || (!dist_matrix) || (!Zcpy)) {
		printf("Out of Memory for C, dist_matrix, Zcpy and Z\n ");
		return -2;
	}

	if (check == 1 && (!Ccpy)) {
		printf("Out of Memory for Ccpy\n ");
		return -2;
	}



	// Generate Observations Vector (Z) for testing phase
	PLASMA_ML_GenerateZVec_tile();


	START_TIMING(time_opt);
	//Maximum Likelihood (Optimization Phase using dlib library)
 	/*max_theta_hat=dlib::find_max_single_variable (
        PLASMA_MLE,
       	starting_theta,
	0.00001, //begin
   	10, //end
        1e-5, //epsilon
        1000, //number of iterations
       	1
	);*/
PLASMA_MLE(0.09);


//PLASMA_MLE(0.9);
	STOP_TIMING(time_opt);



	printf("No. of iteration to converage=%d\n",iter_count);
	printf("Total Optimization Time= %6.2f\n",time_opt);
	printf("Max Theta_hat: %6.2f\n",max_theta_hat);

	//Free Memory Allocation
	PLASMA_ML_finalize();

	return 0;
}





//Generate Observations Vector (Z) for testing Maximum Likelihood function
int PLASMA_ML_GenerateZVec_tile() {

	//Initialization
	double flops=0;
	double eps;
	int info_factorization;


	eps = BLAS_dfpinfo(blas_eps);

	if (verbose == 1)
		fprintf(stderr, "Initializing Distance Matrix ...\n");


	//Uniform random generation of distance matrix (For testing Phase)
	GenerateDistanceMatrix(dist_matrix, N);

	if (verbose == 1)
		fprintf(stderr, "Done ...\n");

	if (verbose == 1)
		fprintf(stderr, "Initializing Co-variance Matrix ...\n");
	//Generate co-variance matrix C
	GenerateCovMatrix(C, dist_matrix, N, THETA);

	if (verbose == 1)
		fprintf(stderr, "Done ...\n");

	if (check == 1)
		PLASMA_dlacpy(PlasmaUpperLower, N, N, C, LDC, Ccpy, LDC);

	//For Matlab Test
	//writetofile("matrix.txt", N, N, C, N); //Sameh;

	if (verbose == 1)
		fprintf(stderr, "Cholesky factorization of Sigma (Generation Phase) ...\n");

	//Cholesky factorization
	int success = PLASMA_dpotrf(PlasmaLower, N, C, LDC);

	if (success != PLASMA_SUCCESS) {
		printf("Factorization cannot be performed..\n"
				"The matrix is not positive definite\n\n");
		exit(0);
	}


	flops = flops + FLOPS_DPOTRF(N);
	if (verbose == 1)
		fprintf(stderr, " Done.\n");


	//**************************TESTING DPOTRF
	if (check == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "\n");
		info_factorization = check_factorization(N, Ccpy, C, LDC, PlasmaLower, eps);
		if (info_factorization == 0) {
			printf("***************************************************\n");
			printf(" ---- TESTING DPOTRF (%s)............ PASSED !\n", "LOWER");
			printf("***************************************************\n");
		}
		fprintf(stderr, "\n\n");
	}
	//*********************************************

	if (verbose == 1)
		fprintf(stderr, "Initialization of measurement vector ...");



	//Generate Z Vector based on C
	int * iseed = (int *) malloc(4 * sizeof(int));
	iseed[0] = 371;
	iseed[1] = 371;
	iseed[2] = 371;
	iseed[3] = 371;

	//Uniform random generation of e --  ei~N(0,1)
	LAPACKE_dlarnv(3, iseed, N, Z);

	//Calculate Z=V.e where C=VV-1
	cblas_dtrmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, N, C,
			LDC, Z, 1);

	//Copy Z to Zcpy
	PLASMA_dlacpy(PlasmaUpperLower, N, NRHS, Z, LDZ, Zcpy, LDZ);


	if (verbose == 1)
		fprintf(stderr, "Done Z Vector Generation Phase.\n");


	return 0;
}

double PLASMA_MLE(double theta_val) {

	//Initialization
	double theta_hat,det=1.0,logdet=0.0;
	double eps;
	int info_factorization;
	int i=0;
	double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0;
	double flops = 0.0;

	eps = BLAS_dfpinfo(blas_eps);


	//Generate new co-variance matrix C based on new theta
	GenerateCovMatrix(C, dist_matrix, N, theta_val);


	if (check == 1)
		PLASMA_dlacpy(PlasmaUpperLower, N, N, C, LDC, Ccpy, LDC);

	//re-store old Z
	PLASMA_dlacpy(PlasmaUpperLower, N, NRHS, Zcpy, LDZ, Z, LDZ);

	

	//Calculate Cholesky Factorization (C=LL-1)
	if (verbose == 1)
		fprintf(stderr, "Cholesky factorization of Sigma...");

	START_TIMING(time_facto);
	int success = PLASMA_dpotrf(PlasmaLower, N, C, LDC);
	STOP_TIMING(time_facto);

	if (success != PLASMA_SUCCESS) {
		printf("Factorization cannot be performed..\n"
				"The matrix is not positive definite\n\n");

		exit(0);
	}


	flops = flops + FLOPS_DPOTRF(N);
	if (verbose == 1)
		fprintf(stderr, " Done.\n");


	//**************************TESTING DPOTRF
	if (check == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "\n");
		info_factorization = check_factorization(N, Ccpy, C, LDC, PlasmaLower, eps);
		if (info_factorization == 0) {
			printf("***************************************************\n");
			printf(" ---- TESTING DPOTRF (%s)............ PASSED !\n", "LOWER");
			printf("***************************************************\n");
		}
		fprintf(stderr, "\n\n");
	}
	///***********************************************************


	//Calculate log(|C|) --> log(square(|L|))
	if (verbose == 1)
		fprintf(stderr, "Calculating the log determinant ...");

	START_TIMING(logdet_calculate);

	for (i = 0; i < N; i++)
		det = det * C[i * N + i];

	logdet= det==0? 0:log(det*det);


	STOP_TIMING(logdet_calculate);


	if (verbose == 1)
		fprintf(stderr, " Done.\n");


	if (verbose == 1)
		printf("log-determinant=%f\n\n", logdet);


	//Solving Linear System (L*X=Z)--->inv(L)*Z
	if (verbose == 1)
		fprintf(stderr, "Solving the linear system ...\n");


	START_TIMING(time_solve);

	//Compute triangular solve LC*X = Z
	PLASMA_dtrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, N, NRHS, 1,
			C, LDC, Z, LDZ); //Sameh  NRHS=1  Result is overwritten over Z

	STOP_TIMING(time_solve);

	flops = flops + FLOPS_DTRSM(PlasmaLeft,N, NRHS);

	if (verbose == 1)
		fprintf(stderr, "Calculating the log likelihood ...");
	theta_hat = -0.5 * cblas_ddot(N, Z, 1, Z, 1) - 0.5 * logdet
			- (double) (N / 2) * log(2 * PI);

	

	if (verbose == 1)
		fprintf(stderr, " Done.\n");


	// Print Iteration Summary
	printf("***************************************************\n");
	printf("logdet: %f, expr1: %f, expr2:%f\n",logdet,(-0.5 * cblas_ddot(N, Z, 1, Z, 1)), ((double) (N / 2) * log(2 * PI)));
	printf(" ---- Theta: %2.6f ---- Theta_hat: %2.6f\n", theta_val,	theta_hat);
	printf(" ---- Facto Time: %6.2f\n", time_facto);
	printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
	printf(" ---- dtrsm Time: %6.2f\n", time_solve);
	printf(" ---- Total Time: %6.2f\n", time_facto + logdet_calculate + time_solve);
	printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + logdet_calculate + time_solve));
	printf("***************************************************\n");


	iter_count++;
	return theta_hat;
}

void PLASMA_ML_finalize() {

	//writetofile("theta_hat.txt", iter, 1, theta_hat_plot,		sizeof(theta) / sizeof(double));
	//writetofile("theta.txt", iter, 1, theta, sizeof(theta) / sizeof(double));
	free(C);
	free(Z);
	free(Zcpy);
	if (check == 1) {
		free(Ccpy);

	}
}

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
static int check_factorization(int N, double *A1, double *A2, int LDA, int uplo,
		double eps) {
	double Anorm, Rnorm;
	double alpha;
	int info_factorization;
	int i, j;

	double *Residual = (double *) malloc(N * N * sizeof(double));
	double *L1 = (double *) malloc(N * N * sizeof(double));
	double *L2 = (double *) malloc(N * N * sizeof(double));
	double *work = (double *) malloc(N * sizeof(double));

	memset((void*) L1, 0, N * N * sizeof(double));
	memset((void*) L2, 0, N * N * sizeof(double));

	alpha = 1.0;

	LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, ' ', N, N, A1, LDA, Residual, N);

	/* Dealing with L'L or U'U  */
	if (uplo == PlasmaUpper) {
		LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'u', N, N, A2, LDA, L1, N);
		LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'u', N, N, A2, LDA, L2, N);
		cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
				N, N, (alpha), L1, N, L2, N);
	} else {
		LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'l', N, N, A2, LDA, L1, N);
		LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'l', N, N, A2, LDA, L2, N);
		cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
				N, N, (alpha), L1, N, L2, N);
	}

	/* Compute the Residual || A -L'L|| */
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			Residual[j * N + i] = L2[j * N + i] - Residual[j * N + i];

	BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, N, Residual, N, &Rnorm);
	BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, N, A1, LDA, &Anorm);

	printf("========================================================\n");
	printf("Checking the Cholesky Factorization \n");
	printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n", Rnorm / (Anorm * N * eps));

	if ( isnan(Rnorm / (Anorm * N * eps)) || isinf(Rnorm / (Anorm * N * eps))
			|| (Rnorm / (Anorm * N * eps) > 60.0)) {
		printf("-- Factorization is suspicious ! \n");
		info_factorization = 1;
	} else {
		printf("-- Factorization is CORRECT ! \n");
		info_factorization = 0;
	}

	free(Residual);
	free(L1);
	free(L2);
	free(work);

	return info_factorization;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */
/*
static int check_solution(int N, int NRHS, double *A1, int LDA, double *B1,
		double *B2, int LDB, double eps) {
	int info_solution;
	double Rnorm, Anorm, Xnorm, Bnorm;
	double alpha, beta;
	double *work = (double *) malloc(N * sizeof(double));

	eps = LAPACKE_dlamch_work('e');

	alpha = 1.0;
	beta = -1.0;

	Xnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N,
			NRHS, B2, LDB, work);
	Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N,
			A1, LDA, work);
	Bnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N,
			NRHS, B1, LDB, work);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, (alpha), A1,
			LDA, B2, LDB, (beta), B1, LDB);
	Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N,
			NRHS, B1, LDB, work);

	printf("============\n");
	printf("Checking the Residual of the solution \n");
	printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n",
			Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps));

	if ( isnan(Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps))
			|| (Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps) > 10.0)) {
		printf("-- The solution is suspicious ! \n");
		info_solution = 1;
	} else {
		printf("-- The solution is CORRECT ! \n");
		info_solution = 0;
	}

	free(work);

	return info_solution;
}
*/
//Sameh
/* Auxiliary routine: printing a matrix */
/*
void print_matrix(char* desc, int m, int n, double* a, int lda) {
	int i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			printf(" %6.4f", a[j * lda + i]);
		printf("\n");
	}
}

void writetofile(char* fname, int m, int n, double* a, int lda) {
	int i, j;
	FILE *file = fopen(fname, "w");
	for (i = 0; i < m; i++) {
		for (j = 0; j < n - 1; j++)
			fprintf(file, "%6.2f,", a[j * lda + i]);
		fprintf(file, "%6.2f", a[j * lda + i]);
		fprintf(file, "\n");
	}

	fclose(file);
}

*/

double uniform_distribution(double rangeLow, double rangeHigh) {
	double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
	double range = rangeHigh - rangeLow;
	double myRand_scaled = (myRand * range) + rangeLow;
	return myRand_scaled;
}

double CalculateDistance(double x1, double y1, double x2, double y2) {
	double distance = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
	return distance;
}

int GenerateDistanceMatrix(double* dist_matrix, int n) {



	double * X = (double *) malloc(n * sizeof(double));
	double * Y = (double *) malloc(n * sizeof(double));


	int i, j;

	for (i = 0; i < n; i++) {
		X[i] = (R - 0.5 + uniform_distribution(MIN_RAND, MAX_RAND))
					* sqrt((double) n);
		Y[i] = (L - 0.5 + uniform_distribution(MIN_RAND, MAX_RAND))
					* sqrt((double) n);
	}


	printf("n=%d\n",n);

	for (i = 0; i < n; i++) {

		dist_matrix[i * n + i] = 1;  //theta_1=1
		for (j = i + 1; j < n; j++) {
			double dist = CalculateDistance(X[i], Y[i], X[j], Y[j]);
			dist_matrix[i * n + j] = dist;
			dist_matrix[j * n + i] = dist;

		}
	}

	return 1;
}

int GenerateCovMatrix(double* C, double* dist_matrix, int N, double theta) {
	int i;
	for (i = 0; i < N * N; i++) {

		C[i] = exp(-dist_matrix[i] / theta);

	}

	for (i = 0; i < N; i++) {
		C[i * N + i] = 1;  //theta_1=1
	}

	return 1;
}

