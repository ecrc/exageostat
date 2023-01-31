/**
 *
 * @file testing_dposv.c
 *
 *  LAPACKE testing routines
 *  LAPACKE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Sameh Abdulah
 * @date 2021-01-24
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <gsl/gsl_sf_bessel.h>
#include </home/abdullsm/develop/hicma/chameleon/build/install_dir/include/coreblas/cblas.h>
#include <lapacke.h>
#include </home/abdullsm/develop/hicma/chameleon/build/install_dir/include/coreblas/coreblas.h>
#include "flops.h"
#include <nlopt.h>
#include "plasma_types.h"
//#include<mkl.h>
//Sameh Testing MOde
#define MIN_RAND       -0.4
#define MAX_RAND        0.4
#define R        2
#define L       3
#define USAGE(args, details)                       \
  printf( " Proper Usage is : ./main_genton "args" \n" \
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

static int GenerateDistanceMatrix(int n);
static int GenerateCovMatrix(double* C, double * X, double * Y, int N, double *theta);
static double CalculateDistance(double x1, double y1, double x2, double y2);
static double uniform_distribution(double rangeLow, double rangeHigh);
//static void writetofile(char* fname, int m, int n, double* a, int lda);
//static void print_matrix(char* desc, int m, int n, double* a, int lda);

int LAPACKE_ML_GenerateZVec_tile(double * initial_theta);
double LAPACKE_MLE(unsigned n, const double * theta, double * grad, void * my_func_data);
void LAPACKE_ML_finalize();
void write_to_file(char * path, int matrix_size,int nthreads, double gflops, double avg_time);

double * C ;
double * XLOCs,*YLOCs;
double * Z;
double * Zcpy;
double * Ccpy;
size_t N,NRHS,LDC,LDZ,check,verbose;
int iter_count = 0;
int ncores;
int num_params = 6;
int main(int argc, char **argv) {

    /* Check for number of arguments*/
    if (argc != 8) {
        USAGE("NCORES N TS CHECK VERBOSE",
                "   - NCORES   : number of cores\n" "   - N        : the size of the matrix\n" //sameh
                "   - CHECK    : check the factorization and the solution\n" "   - VERBOSE  : verbose\n");
        return -1;
    }

    ncores = atoi(argv[1]);
    N = atoi(argv[2]);
    //int ts = atoi(argv[3]);
    check = atoi(argv[3]);
    verbose = atoi(argv[4]);
    double time_opt=0.0;
    double max_theta_hat=0.0;
    NRHS = 1;
    LDC = N;
    LDZ = N;

    printf("================================================================================================ N:%d\n",N);

    //Memory Allocation
    C = (double *) malloc((size_t)LDC *(size_t) N * sizeof(double));
    Z = (double *) malloc(LDZ * (size_t)NRHS * sizeof(double));
    Zcpy = (double *) malloc(LDZ * (size_t)NRHS * sizeof(double));
    XLOCs = (double *) malloc((size_t)N * sizeof(double));
    YLOCs = (double *) malloc((size_t)N * sizeof(double));

    if (check == 1)
        Ccpy = (double *) malloc(LDC * N * sizeof(double));

    /* Check if unable to allocate memory */
    if ((!C) || (!Z) ||  (!Zcpy)) {
        printf("Out of Memory for C,  Zcpy and Z\n ");
        return -2;
    }

    if (check == 1 && (!Ccpy)) {
        printf("Out of Memory for Ccpy\n ");
        return -2;
    }

    //stop gsl error handler
    gsl_set_error_handler_off () ;

    double * initial_theta=(double *) malloc(num_params * sizeof(double));
    
    initial_theta[0]=atof(argv[5]);
    initial_theta[1]=atof(argv[6]);
    initial_theta[2]=atof(argv[7]);
    initial_theta[3]=atof(argv[8]);
    initial_theta[4]=atof(argv[9]);
    initial_theta[5]=atof(argv[10]);


    // Generate Observations Vector (Z) for testing phase
    LAPACKE_ML_GenerateZVec_tile(initial_theta);


    double * starting_theta=(double *) malloc(num_params * sizeof(double));
    starting_theta[0]=1;
    starting_theta[1]=0.1;
    starting_theta[2]=0.5;
    starting_theta[3]=1;
    starting_theta[4]=0.1;
    starting_theta[5]=0.5;    

    nlopt_opt opt;
    opt=nlopt_create( NLOPT_LN_BOBYQA, num_params);
    double lb[6]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    double up[6]={5, 5, 5, 4, 5, 5};
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, up);
    double * opt_f;
    nlopt_set_xtol_rel(opt, 1e-5);	


    START_TIMING(time_opt);

    nlopt_set_max_objective(opt, LAPACKE_MLE ,NULL );
    nlopt_optimize(opt, starting_theta, &opt_f);
    //Maximum Likelihood (Optimization Phase using dlib library)
    //LAPACKE_MLE(0.5);
    //LAPACKE_MLE(0.9);
    STOP_TIMING(time_opt);



    printf("No. of iteration to converage=%d\n",iter_count);
    printf("Total Optimization Time= %6.2f\n",time_opt);
    printf("Max Theta_hat: %6.2f\n",max_theta_hat);

    //Free Memory Allocation
    LAPACKE_ML_finalize();

    return 0;
}





//Generate Observations Vector (Z) for testing Maximum Likelihood function
int LAPACKE_ML_GenerateZVec_tile(double * initial_theta) {

    //Initialization
    double flops=0;
    double eps;
    int info_factorization;
    double time_gen=0.0;
    double time_dportf=0.0;


    eps = BLAS_dfpinfo(blas_eps);

    if (verbose == 1)
        fprintf(stderr, "Initializing Distance Matrix ...\n");


    //Uniform random generation of distance matrix (For testing Phase)
    GenerateDistanceMatrix(N);

    //printf(" ---- Matrix Generation Time: %6.2f\n", time_gen);


    if (verbose == 1)
        fprintf(stderr, "Done ...\n");

    START_TIMING(time_gen);
    if (verbose == 1)
        fprintf(stderr, "Initializing Co-variance Matrix ...\n");
    //Generate co-variance matrix C
    GenerateCovMatrix(C, XLOCs,YLOCs, N, initial_theta);
    STOP_TIMING(time_gen);
    /*
       int index2=0;
       for (index2=0;index2 <16*16 ;index2++)
       {        	printf("%f ",C[index2],C[index2]);
       if((index2+1)%16==0)
       printf("\n");
       }
       printf("\n");
       exit(0);
       */
    printf(" ---- Matrix Generation Time: %6.2f\n", time_gen);


    if (verbose == 1)
        fprintf(stderr, "Done ...\n");



    if (check == 1)
        LAPACKE_dlacpy(LAPACK_COL_MAJOR,'L', N, N, C, LDC, Ccpy, LDC);

    //For Matlab Test
    //writetofile("matrix.txt", N, N, C, N); //Sameh;

    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma (Generation Phase) ...\n");

    //Cholesky factorization
    START_TIMING(time_dportf);

    int success = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L', N, C, LDC);

    STOP_TIMING(time_dportf);

    printf("dportf execution time: %f\n",time_dportf);
    /*	if (!success) {
        printf("Factorization cannot be performed..\n"
        "The matrix is not positive definite\n\n");
        exit(0);
        }
        */

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
    iseed[0] = 0;
    iseed[1] = 0;
    iseed[2] = 0;
    iseed[3] = 1;

    //Uniform random generation of e --  ei~N(0,1)
    LAPACKE_dlarnv(3, iseed, N, Z);

    //Calculate Z=V.e where C=VV-1
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
            N, NRHS, 1, C, LDC, Z, LDZ);

    //Copy Z to Zcpy
    LAPACKE_dlacpy(LAPACK_COL_MAJOR,'L', N, NRHS, Z, LDZ, Zcpy, LDZ);

    //int index=0;
    //for (index=0;index <16;index++)
    //printf("%d - ",Z[index]);
    //exit(0);

    if (verbose == 1)
        fprintf(stderr, "\nDone Z Vector Generation Phase.\n");


    return 0;
}

double LAPACKE_MLE(unsigned n, const double * theta, double * grad, void * my_func_data) {

    //Initialization
    double theta_hat,det=1.0,logdet=0.0;
    double eps;
    int info_factorization;
    int i=0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, time_gen=0.0;
    double flops = 0.0;

    eps = BLAS_dfpinfo(blas_eps);

    START_TIMING(time_gen);

    //Generate new co-variance matrix C based on new theta
    GenerateCovMatrix(C, XLOCs, YLOCs, N, (double*)theta);

    STOP_TIMING(time_gen);



    printf("\n Matrix Generation Time: %6.2f\n", time_gen);

    if (check == 1)
        LAPACKE_dlacpy(CblasColMajor,PlasmaGeneral, N, N, C, LDC, Ccpy, LDC);

    //re-store old Z
    LAPACKE_dlacpy(CblasColMajor,PlasmaGeneral, N, NRHS, Zcpy, LDZ, Z, LDZ);



    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "--Cholesky factorization of Sigma...");
    //printf ("N: %d\n",N);

    //int k=0;
    //	for (k=0;k<N;k++)
    //		printf("%f - ", C[k]);
    //	exit(0);


    START_TIMING(time_facto);
    int success = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L', N, C, LDC);
    STOP_TIMING(time_facto);
    /*
       if (!success ) {
       printf("Factorization cannot be performed..\n"
       "The matrix is not positive definite\n\n");

       exit(0);
       }
       */

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

    logdet= det==0?0:2*log(det);

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
    cblas_dtrsm(CblasColMajor,CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, N, NRHS, 1, C, LDC, Z, LDZ); //Sameh  NRHS=1  Result is overwritten over Z

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
    printf("logdet: %f, dotp: %f, expr2:%f\n",logdet,cblas_ddot(N, Z, 1, Z, 1), ((double) (N / 2) * log(2 * PI)));
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Total Time: %6.2f\n", time_facto + logdet_calculate + time_solve+time_gen);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + time_solve));
    printf("***************************************************\n");


    write_to_file("results.txt",N,64, (flops / 1e9 / (time_facto + time_solve)), (time_facto + logdet_calculate + time_solve+time_gen));

    exit(0);
    iter_count++;
    return theta_hat;
}

void LAPACKE_ML_finalize() {

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
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
    //return distance;
}

int GenerateDistanceMatrix( int n) {
    int r=0, l=0;
    int index=0;
    int sqrtn=sqrt(n);
    double * Xt = (double *) calloc((size_t)sqrtn, sizeof(double));
    double * Yt = (double *) calloc((size_t)sqrtn , sizeof(double));

    for(r=0;r<sqrtn;r++)
        Xt[r] = ((r+1) - 0.5 + uniform_distribution(-0.4, 0.4))/sqrtn;

    for(l=0;l<sqrtn;l++)
        Yt[l] = ((l+1) - 0.5 + uniform_distribution(-0.4, 0.4))/sqrtn;


    for(r=0;r<sqrtn;r++)
        for(l=0;l<sqrtn;l++){
            XLOCs[index]=Xt[r];
            YLOCs[index++]=Yt[l];
        }

    free(Xt);
    free(Yt);	
    return 1;
}

static int GenerateCovMatrix(double* C, double* X, double * Y, int N, double* theta) {
    size_t  i;
    //double value;
    //double x0,y0;
    //double expr;
    double con=0.0;

    con=pow(2,(theta[2]-1)) * tgamma(theta[2]);
    con=1.0/con;
    con=theta[0]*con;	

    //#pragma omp parallel for private(i, j, expr) schedule(dynamic) num_threads(ncores) 
#pragma omp parallel for schedule(dynamic) num_threads(ncores) 
    for (i = 0; i < N; i++) {
        double x0=XLOCs[i];
        double y0=YLOCs[i];

        size_t j;

        for (j = 0 ; j < N; j++) {

            double expr=CalculateDistance(x0, y0, XLOCs[j], YLOCs[j] )/theta[1];
            if(expr==0)
                C[j * N + i]=theta[0];
            else
                C[j * N + i]=con*pow(expr,theta[2])*gsl_sf_bessel_Knu(theta[2],expr); // Matern Function

            //	C[i * N + j]=C[j * N + i];
        }
    }

    //exit(0);
    return 1;
}


void write_to_file(char * path, int matrix_size,int nthreads, double gflops, double avg_time)
{

    FILE *pFile;

    pFile=fopen(path,"a");

    if(pFile==NULL) {

        printf("Cannot access the results path\n");
        exit(0);
    }

    fprintf(pFile, "%d\t", matrix_size);
    fprintf(pFile, "%d\t", nthreads);
    fprintf(pFile, "%f\t", gflops);
    fprintf(pFile, "%f\n", avg_time);
    fclose(pFile);
}


