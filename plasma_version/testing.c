/**
 *
 * @file testing.c
 *
 *  
 *  MLE is a software package provided by KAUST,
 *  King Abdullah University of Science and Technology - ECRC
 *
 * @version 1.0.0
 * @author Sameh Abdulsh
 * @date 2016-11-22
 * @generated d Fri Nov 22 15:11:13 2016
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
//#include <dlib/optimization.h>
#include <nlopt.h>
#include <math.h>
#include <plasma.h>
#include <morse.h>
#include <starpu.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "flops.h"
#include <argp.h>

const char *argp_program_version = "Version 1.0";
const char *argp_program_bug_address = "<sameh.abdulah@kaust.edu.sa>";
struct arguments {
//  char *args[0];          /* ncores, N, ts, locs_file, obs_dir, and obs_timestamp */
    int test;                   /* The -t flag (USE Synthetic Dataset*/
    int check;                /* The -c flag */
    int verbose;             /* The -v flag */
    char *ncores;              //Number of cores
    char *N;                   //Matrix size (Only in the case oif sysnthetic dataset)
    char *ts;                  //Tiling size
    char *kernel;
    int chameleon;          // 0--> plasma 1-->chameleon  The -sys flag
    int async;          // 0-->tile  1-->tile_async       The -async flag
    char *locs_file;        // Locations files (in the case of real dataset
    char *obs_dir;          // Observations directory in the case of real dataset
    int timestamp;       // obstervations timestamp

};

static struct argp_option options[] =
        {
                {"test",      't', 0,                        0, "Execute in test mode"},
                {"check",     'c', 0,                        0, "Produce check output"},
                {"verbose",   'v', 0,                        0, "Produce verbose output"},
                {"ncores",    'n', "NCORES",                 0, "Number of Cores"},
                {"N",         's', "MATRIX_SIZE",            0, "Synthetic Matrix Size"},
                {"ts",        'e', "TILE_SIZE",              0, "Number of Tiles"},
                {"kernel",    'k', "KERNEL",                 0, "Computation Model"},
                {"chameleon", 'b', 0,                        0, "Use Chameleon Instead of Plasma"},
                {"async",     'a', 0,                        0, "Asynchronous"},
                {"locs_file", 'l', "LOCATIONS_FILE",         0, "Read Locations from this Location File"},
                {"obs_dir",   'o', "OBSERVATIONS_DIRECTORY", 0, "Read Observations from this Directory Path"},
                {"timestamp", 'p', "TIMESTAMP",              0, "Observation Timestamp"},
                {0}
        };


static error_t
parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;

    switch (key) {
        case 't':
            arguments->test = 1;
            break;
        case 'c':
            arguments->check = 1;
            break;
        case 'v':
            arguments->verbose = 1;
            break;
        case 'k':
            arguments->kernel = arg;
            break;

        case 'n':
            arguments->ncores = arg;  //non-optional;
            break;
        case 's':
            arguments->N = arg;
            break;
        case 'e':
            arguments->ts = arg;  //non-optional
            break;
        case 'b':
            arguments->chameleon = 1;
            break;
        case 'a':
            arguments->async = 1;
            break;
        case 'l':
            arguments->locs_file = arg;
            break;
        case 'o':
            arguments->obs_dir = arg;
            break;
        case 'p':
            arguments->timestamp = arg;
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static char args_doc[] = "";

static char doc[] =
        "NLE -- A program to evaluate Maximum Likelihood function using both real and synthetic dataset on two different platforms (PLASMA (default) and Chameleon)";

static struct argp argp = {options, parse_opt, args_doc, doc};


// For Bessel Function (Do not change)
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
//tiling
#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )
#define mBLKLDD(A, k) A->get_blkldd( A,k )
#define A(m, n) (double *)plasma_getaddr(A, m, n)
//Sameh Testing MOde
#define MIN_RAND       -0.4
#define MAX_RAND        0.4
#define R        2
#define L       3
#define THETA1 0.1 //weak=0.03, medium=0.1,  strong=0.3
#define THETA2 0.03
#define THETA3 0.03
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


static double CalculateDistance(double x1, double y1, double x2, double y2);

static double uniform_distribution(double rangeLow, double rangeHigh);

//static void writetofile(char* fname, int m, int n, double* a, int lda);
static void print_matrix(char *desc, int m, int n, double *a, int lda);

int PLASMA_MLE_GenZVec_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
                                  PLASMA_request *request);

int MORSE_MLE_GenZVec_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence,
                                 MORSE_request_t *request);


static int GenerateXYLoc(int n, char *locs_file);


static double CalculateDistance(double x1, double y1, double x2, double y2);

static double uniform_distribution(double rangeLow, double rangeHigh);

//tiling
static void CORE_dGenCovMat_quark(Quark *quark);

static void CORE_dGenZVec_quark(Quark *quark);

static int PLASMA_MLE_GenCovMat_Tile_Async(PLASMA_desc *descA,
                                           PLASMA_sequence *sequence, PLASMA_request *request, double *theta);

static int MORSE_MLE_GenCovMat_Tile_Async(MORSE_desc_t *descA,
                                          MORSE_sequence_t *sequence, MORSE_request_t *request, double *theta);


static void QUARK_CORE_dGenCovMat(Quark *quark, Quark_Task_Flags *task_flags,
                                  int m, int n, double *A, int lda, int m0, int n0, double *theta);

static void matcov_comp_Tile(double *A, int m, int n, int m0, int n0, double *theta);

static void CORE_det_quark(Quark *quark);

static void det_comp_Tile(double *A, int m, int n, int m0, int n0, double *determinant);


static void QUARK_CORE_det(Quark *quark, Quark_Task_Flags *task_flags, int m,
                           int n, double *A, int lda, int m0, int n0, double *determinant);

static int PLASMA_MLE_det_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
                                     PLASMA_request *request, double *result);

static int MORSE_MLE_det_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence,
                                    MORSE_request_t *request, double *result);

void QUARK_CORE_dGenZVec(Quark *quark, Quark_Task_Flags *task_flags, int m,
                         int n, double *A, int lda, int m0, int n0, double *r);

void veccov_comp_Tile(double *Z, int m, int n, int m0, int n0, double *r);

int PLASMA_ML_GenerateZVec_Tile(char *obs_di, int obs_timestamp);


double PLASMA_MLE_Tile(unsigned n, const double *theta, double *grad, void *my_func_data);


void ML_finalize();

int PLASMA_ML_GenerateZVec_Tile_Async(char *obs_dir, int obs_timestamp);

double PLASMA_MLE_Tile_Async(unsigned n, const double *theta, double *grad, void *my_func_data);

int MORSE_ML_GenerateZVec_Tile(char *obs_dir, int obs_timestamp);

double MORSE_MLE_Tile(unsigned n, const double *theta, double *grad, void *my_func_data);

int MORSE_ML_GenerateZVec_Tile_Async(char *obs_dir, int obs_timestamp);


double MORSE_MLE_Tile_Async(unsigned n, const double *theta, double *grad, void *my_func_data);

int countlines(char *filename);

static double bessi0(double x);

static double bessi1(double x);


double bessk(int n, double x);

static double bessk0(double x);

static double bessk1(double x);


double *C;
double *Z;
double *Zcpy;
double *Ccpy;
//********************************************************************
PLASMA_desc *descC = NULL;
PLASMA_desc *descZ = NULL;
PLASMA_desc *descZcpy = NULL;
MORSE_desc_t *MORSE_descC = NULL;
MORSE_desc_t *MORSE_descZ = NULL;
MORSE_desc_t *MORSE_descZcpy = NULL;
double *X;
double *Y;
//PLASMA sequence uniquely identifies a set of asynchronous function calls sharing common exception handling.
PLASMA_sequence *sequence;
//PLASMA request uniquely identifies each asynchronous function call.
PLASMA_request request[19] = {PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,

                              PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER};
//MORSE sequence uniquely identifies a set of asynchronous function calls sharing common exception handling.
MORSE_sequence_t *msequence;
//MORSE request uniquely identifies each asynchronous function call.
MORSE_request_t mrequest[19] = {MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER, MORSE_REQUEST_INITIALIZER};
//***************************************************************************************

int N, NRHS, LDC, LDZ, async, check, verbose;
int iter_count = 0;
int based_sys;
double *THETA;
char *kernel;

int main(int argc, char **argv) {

    /* Check for number of arguments*/
/*	if (argc != 11) {
		USAGE("NCORES N TS CHECK VERBOSE",
				"   - NCORES   : number of cores\n" "   - N        : the size of the matrix\n" "   - TS       : tile size\n" //sameh
				"   - BASED_SYS : Plasma (0) Chameleon(1) \n- ASYNC: asynchronous \n  - CHECK    : check the factorization and the solution\n" "   - VERBOSE  : verbose  \n");
		return -1;
	}

	int ncores = atoi(argv[1]);
	N = atoi(argv[2]);
	int ts = atoi(argv[3]);
        based_sys = atoi(argv[4]); // 0--> plasma 1-->chameleon
	async= atoi(argv[5]); // 0-->tile  1-->tile_async
        char * locs_file= argv[6];   
 	char * obs_dir= argv[7];
	int obs_timestamp=atoi(argv[8]);
	check = atoi(argv[9]);
	verbose = atoi(argv[10]);

*/


// new

    struct arguments arguments;
    FILE *outstream;
    arguments.test = 0;                  /* execute in test mode */
    arguments.check = 0;                 /* The -c flag */
    arguments.verbose = 0;             /* The -v flag */
    arguments.chameleon = 0;           // 0--> plasma 1-->chameleon  The -sys flag
    arguments.async = 0;               // 0-->tile  1-->tile_async       The -async flag
    arguments.kernel = "EXP_COV_MODEL"; // EXP_COV_MODEL or MATERN_MODEL
    arguments.ncores = "1";              //Number of cores
    arguments.N = "0";                   //Matrix size (Only in the case oif sysnthetic dataset)
    arguments.ts = "0";                  //Tiling size
    arguments.locs_file = "";        // Locations files (in the case of real dataset
    arguments.obs_dir = "";          // Observations directory in the case of real dataset
    arguments.timestamp = 0;       // obstervations timestamp



    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    int test = arguments.test;
    int ncores = atoi(arguments.ncores);
    N = atoi(arguments.N);
    kernel = arguments.kernel;
    int ts = atoi(arguments.ts);
    based_sys = arguments.chameleon;// 0--> plasma 1-->chameleon
    async = arguments.async; // 0-->tile  1-->tile_async
    char *locs_file = arguments.locs_file;
    char *obs_dir = arguments.obs_dir;
    int obs_timestamp = arguments.timestamp;
    check = arguments.check;
    verbose = arguments.verbose;


    if (test == 0 && (strcmp(locs_file, "") == 0 || strcmp(obs_dir, "") == 0 || obs_timestamp == 0 || N > 0)) {
        printf("\nIn real mode: please use locs_file, obs_dir, and obs_timestamp arguments only and ignore N (Matrix Size)\n\n");
        exit(0);
    } else if (test == 1 && (strcmp(locs_file, "") != 0 || strcmp(obs_dir, "") != 0 || obs_timestamp != 0 || N == 0)) {
        printf("\nIn test mode: please use N (Matrix Size) and ignore locs_file, obs_dir, and obs_timestamp arguments\n\n");
        exit(0);

    }

    if (strcmp(kernel, "EXP_COV_MODEL") != 0 && strcmp(kernel, "MATERN_MODEL") != 0) {
        printf("\nPlease use EXP_COV_MODEL or MATERN_MODEl for kernel input... default value is (EXP_COV_MODEL)\n\n");
        exit(0);
    }



//*****************************************************************************************************************************

    int locs;

    if (strcmp(locs_file, "") != 0) {
        printf("locs file: %s\n", locs_file);
        locs = countlines(locs_file);
        N = locs;
    }

    double time_opt = 0.0;


    //matern function
    double *starting_theta = (double *) malloc(3 * sizeof(double));
    starting_theta[0] = 0.03;
    starting_theta[1] = 0.03;
    starting_theta[2] = 0.03;
    THETA = (double *) malloc(3 * sizeof(double));
    THETA[0] = THETA1;
    THETA[1] = THETA2;
    THETA[2] = THETA3;


    double max_theta_hat = 0.0;
    NRHS = 1;
    LDC = N;
    LDZ = N;

    //Memory Allocation
    X = (double *) malloc(N * sizeof(double));
    Y = (double *) malloc(N * sizeof(double));
    C = (double *) malloc(LDC * N * sizeof(double));
    Z = (double *) malloc(LDZ * NRHS * sizeof(double));
    Zcpy = (double *) malloc(LDZ * NRHS * sizeof(double));



    /* Check if unable to allocate memory */
    if ((!C) || (!Z) || (!Zcpy) || (!X) || (!Y)) {
        printf("Out of Memory for C,  Zcpy , Z, X, and Y\n ");
        return -2;
    }


    GenerateXYLoc(N, locs_file);


    nlopt_opt opt;
    //Initial nlopt (Optimization)
    if (strcmp(kernel, "MATERN_MODEL") == 0)
        opt = nlopt_create(NLOPT_LN_BOBYQA, 3);
    else
        opt = nlopt_create(NLOPT_LN_BOBYQA, 1);

    double lb[3] = {0.001, 0.001, 0.001};
    double up[3] = {10, 10, 10};
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, up);
    double *opt_f;
    nlopt_set_xtol_rel(opt, 1e-5);


    if (based_sys == 0) {
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


        //Identifies a set of routines sharing common exception handling.
        PLASMA_Sequence_Create(&sequence);

        //Create matrix descriptor -- ts -->number of rows in tile , ts -->number of columns in tile, ts*ts-->title size (number of elements), N -->number of rows in the whole matrix, N --> number pf columns in the whole matrix
        PLASMA_Desc_Create(&descC, C, PlasmaRealDouble, ts, ts, ts * ts, LDC, N, 0,
                           0, LDC, N);

        PLASMA_Desc_Create(&descZ, Z, PlasmaRealDouble, ts, ts, ts * ts, LDZ, NRHS, 0, 0, LDZ, NRHS);

        PLASMA_Desc_Create(&descZcpy, Zcpy, PlasmaRealDouble, ts, ts, ts * ts, LDZ, NRHS, 0, 0, LDZ, NRHS);


        // Generate Observations Vector (Z) for testing phase
        if (async == 0) {
            PLASMA_ML_GenerateZVec_Tile(obs_dir, obs_timestamp);
            START_TIMING(time_opt);

            //Maximum Likelihood (Optimization Phase using nlopt library)
            nlopt_set_max_objective(opt, PLASMA_MLE_Tile, NULL);
            nlopt_optimize(opt, starting_theta, &opt_f);
            //printf("success: %d\n", success);
            /*		max_theta_hat=dlib::find_max_single_variable (
                        PLASMA_MLE_Tile,
                        starting_theta,
                        0.00001, //begin
                        10, //end
                        1e-5, //epsilon
                        1000, //number of iterations
                        1
                    );
        */
//			double x=0.09;
            //                      PLASMA_MLE_Tile(0,&x,NULL,NULL);

            STOP_TIMING(time_opt);
        } else {
            PLASMA_ML_GenerateZVec_Tile_Async(obs_dir, obs_timestamp);

            START_TIMING(time_opt);
            nlopt_set_max_objective(opt, PLASMA_MLE_Tile_Async, NULL);
            nlopt_optimize(opt, starting_theta, &opt_f);
            /* //Maximum Likelihood (Optimization Phase using dlib library)
            max_theta_hat=dlib::find_max_single_variable (
                PLASMA_MLE_Tile_Async,
                starting_theta,
                0.00001, //begin
                10, //end
                1e-5, //epsilon
                1000, //number of iterations
                1
            );
            */
//			PLASMA_MLE_Tile_Async(0.9);

            STOP_TIMING(time_opt);
        }

    }

// MORSE based- system
    else {

        /*-------------------------------------------------------------
             *              Initialization of Morse
             */
        MORSE_Init(ncores, 0);  //MORSE_init(NCPU,NGPU);
        MORSE_Disable(PLASMA_AUTOTUNING);
        MORSE_Set(PLASMA_TILE_SIZE, ts);
        MORSE_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING);
        MORSE_Enable(PLASMA_WARNINGS);
        MORSE_Enable(PLASMA_ERRORS);
        /*-------------------------------------------------------------*/


        //Identifies a set of routines sharing common exception handling.
        MORSE_Sequence_Create(&msequence);

        //Create matrix descriptor -- ts -->number of rows in tile , ts -->number of columns in tile, ts*ts-->title size (number of elements), N -->number of rows in the whole matrix, N --> number pf columns in the whole matrix
        MORSE_Desc_Create(&MORSE_descC, C, MorseRealDouble, ts, ts, ts * ts, LDC, N, 0,
                          0, LDC, N, 1, 1);

        MORSE_Desc_Create(&MORSE_descZ, Z, MorseRealDouble, ts, ts, ts * ts, LDZ, NRHS, 0, 0, LDZ, NRHS, 1, 1);

        MORSE_Desc_Create(&MORSE_descZcpy, Zcpy, MorseRealDouble, ts, ts, ts * ts, LDZ, NRHS, 0, 0, LDZ, NRHS, 1, 1);



        // Generate Observations Vector (Z) for testing phase
        if (async == 0) {
            MORSE_ML_GenerateZVec_Tile(obs_dir, obs_timestamp);
            START_TIMING(time_opt);
            nlopt_set_max_objective(opt, MORSE_MLE_Tile, NULL);
            nlopt_optimize(opt, starting_theta, &opt_f);
            //Maximum Likelihood (Optimization Phase using dlib library)
            /* max_theta_hat=dlib::find_max_single_variable (
                            MORSE_MLE_Tile,
                            starting_theta,
                            0.00001, //begin
                            10, //end
                            1e-5, //epsilon
                            1000, //number of iterations
                            1
                );
            */
            //MORSE_MLE_Tile(0.9);
            STOP_TIMING(time_opt);
        } else {
            MORSE_ML_GenerateZVec_Tile_Async(obs_dir, obs_timestamp);
            START_TIMING(time_opt);
            nlopt_set_max_objective(opt, MORSE_MLE_Tile_Async, NULL);
            nlopt_optimize(opt, starting_theta, &opt_f);
            //Maximum Likelihood (Optimization Phase using dlib library)
            /* max_theta_hat=dlib::find_max_single_variable (
                        MORSE_MLE_Tile_Async,
                        starting_theta,
                        0.00001, //begin
                        10, //end
                        1e-5, //epsilon
                        1000, //number of iterations
                        1
            );
    */
            //	MORSE_MLE_Tile_Async(0.9);
            STOP_TIMING(time_opt);
        }


    }

    printf("No. of iteration to converage=%d\n", iter_count);
    printf("Total Optimization Time= %6.2f\n", time_opt);
    printf("Max Theta_hat: %6.2f\n", max_theta_hat);

    //Free Memory Allocation
    ML_finalize();
    nlopt_destroy(opt);

    return 0;
}


//Generate Observations Vector (Z) for testing Maximum Likelihood function
int PLASMA_ML_GenerateZVec_Tile(char *obs_dir, int obs_timestamp) {

    //Initialization
    double flops = 0;


    //Set C diagonal equals to 1 (simple case)
    if (async == 0)

        PLASMA_dlaset_Tile(PlasmaUpperLower, 0, 1, descC);


    if (verbose == 1)
        fprintf(stderr, "Initializing Co-variance Matrix ...\n");


    //Generate co-variance matrix C
    PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], THETA);

    if (verbose == 1)
        fprintf(stderr, "Done ...\n");




    //    int IPV;
    //  PLASMA_dgetrf_Tile(descC,&IPV);
//	double res=0;
    //      PLASMA_dgecon_Tile(PlasmaOneNorm,descC,PlasmaOneNorm,&res);
//	printf("res: %f\n",res);
//	exit(0);


    if (strcmp(obs_dir, "") == 0) {

        if (verbose == 1)
            fprintf(stderr,
                    "Cholesky factorization of Sigma (Generation Phase) ...\n");


        //Cholesky factorization
        int success = PLASMA_dpotrf_Tile(PlasmaLower, descC);

        if (success != PLASMA_SUCCESS) {
            printf("Factorization cannot be performed..\n"
                   "The matrix is not positive definite\n\n");
            exit(0);
        }

        flops = flops + FLOPS_DPOTRF(N);
        if (verbose == 1)
            fprintf(stderr, " Done.\n");

        // Generate Z Vector
        PLASMA_MLE_GenZVec_Tile_Async(descZ, sequence, &request[1]);

        PLASMA_dtrmm_Tile(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1, descC, descZ);
    } else {

        if (verbose == 1)
            fprintf(stderr,
                    "Reading omebservations from dir ...\n");

        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        char *pch;
        int str_count = 0;
        char *zeroes;
        char integer_string[32];
        int i = 0;

        char *obs_dir_name = (char *) malloc(50 * sizeof(char));

        for (i = 0; i < 103; i++) {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i + 1));

            if (i < 9)
                zeroes = "00";
            else if (i >= 9 && i < 99)
                zeroes = "0";
            else
                zeroes = "";

            strcat(obs_dir_name, "/ppt.complete.Y");

            strcat(obs_dir_name, zeroes);
            strcat(obs_dir_name, integer_string);


            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count = 0;

                pch = strtok(line, " ,");
                while (pch != NULL) {
                    if (str_count == obs_timestamp)
                        Z[i] = atof(pch);
                    pch = strtok(NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }



    //Copy Z to Zcpy
    PLASMA_dlacpy_Tile(PlasmaUpperLower, descZ, descZcpy);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");


    return 0;

}


//Generate Observations Vector (Z) for testing Maximum Likelihood function
int PLASMA_ML_GenerateZVec_Tile_Async(char *obs_dir, int obs_timestamp) {

    //Initialization
    double flops = 0;


    //Set C diagonal equals to 1 (simple case)
    if (async == 0)

        PLASMA_dlaset_Tile_Async(PlasmaUpperLower, 0, 1, descC, sequence, &request[1]);


    if (verbose == 1)
        fprintf(stderr, "Initializing Co-variance Matrix ...\n");


    //Generate co-variance matrix C
    PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], THETA);

    if (verbose == 1)
        fprintf(stderr, "Done ...\n");


    if (strcmp(obs_dir, "") == 0) {

        if (verbose == 1)
            fprintf(stderr,
                    "Cholesky factorization of Sigma (Generation Phase) ...\n");


        //Cholesky factorization
        int success = PLASMA_dpotrf_Tile_Async(PlasmaLower, descC, sequence, &request[1]);

        if (success != PLASMA_SUCCESS) {
            printf("Factorization cannot be performed..\n"
                   "The matrix is not positive definite\n\n");
            exit(0);
        }

        flops = flops + FLOPS_DPOTRF(N);
        if (verbose == 1)
            fprintf(stderr, " Done.\n");

        // Generate Z Vector
        PLASMA_MLE_GenZVec_Tile_Async(descZ, sequence, &request[1]);

        PLASMA_dtrmm_Tile_Async(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1, descC, descZ, sequence,
                                &request[1]);
    } else {

        if (verbose == 1)
            fprintf(stderr,
                    "Reading omebservations from dir ...\n");

        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        char *pch;
        int str_count = 0;
        char *zeroes;
        char integer_string[32];
        int i = 0;

        char *obs_dir_name = (char *) malloc(50 * sizeof(char));

        for (i = 0; i < 103; i++) {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i + 1));

            if (i < 9)
                zeroes = "00";
            else if (i >= 9 && i < 99)
                zeroes = "0";
            else
                zeroes = "";

            strcat(obs_dir_name, "/ppt.complete.Y");

            strcat(obs_dir_name, zeroes);
            strcat(obs_dir_name, integer_string);


            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count = 0;

                pch = strtok(line, " ,");
                while (pch != NULL) {
                    if (str_count == obs_timestamp)
                        Z[i] = atof(pch);
                    pch = strtok(NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }



    //Copy Z to Zcpy
    PLASMA_dlacpy_Tile_Async(PlasmaUpperLower, descZ, descZcpy, sequence, &request[1]);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");


    return 0;


}


//Generate Observations Vector (Z) for testing Maximum Likelihood function
int MORSE_ML_GenerateZVec_Tile(char *obs_dir, int obs_timestamp) {

    //Initialization
    double flops = 0;

    //Set C diagonal equals to 1 (simple case)
    if (async == 0)

        MORSE_dlaset_Tile(MorseUpperLower, 0, 1, MORSE_descC);


    if (verbose == 1)
        fprintf(stderr, "Initializing Co-variance Matrix ...\n");


    //Generate co-variance matrix C
    MORSE_MLE_GenCovMat_Tile_Async(MORSE_descC, msequence, &mrequest[1], THETA);


    if (verbose == 1)
        fprintf(stderr, "Done ...\n");

    if (strcmp(obs_dir, "") == 0) {
        if (verbose == 1)
            fprintf(stderr,
                    "Cholesky factorization of Sigma (Generation Phase) ...\n");

        //Cholesky factorization
        int success = MORSE_dpotrf_Tile(MorseLower, MORSE_descC);

        if (success != PLASMA_SUCCESS) {
            printf("Factorization cannot be performed..\n"
                   "The matrix is not positive definite\n\n");
            exit(0);
        }

        flops = flops + FLOPS_DPOTRF(N);
        if (verbose == 1)
            fprintf(stderr, " Done.\n");


        // Generate Z Vector
        MORSE_MLE_GenZVec_Tile_Async(MORSE_descZ, msequence, &mrequest[1]);

        MORSE_dtrmm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ);
    } else {
        if (verbose == 1)
            fprintf(stderr,
                    "Reading omebservations from dir ...\n");

        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        char *pch;
        int str_count = 0;
        char *zeroes;
        char integer_string[32];
        int i = 0;

        char *obs_dir_name = (char *) malloc(50 * sizeof(char));

        for (i = 0; i < 103; i++) {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i + 1));

            if (i < 9)
                zeroes = "00";
            else if (i >= 9 && i < 99)
                zeroes = "0";
            else
                zeroes = "";

            strcat(obs_dir_name, "/ppt.complete.Y");

            strcat(obs_dir_name, zeroes);
            strcat(obs_dir_name, integer_string);


            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count = 0;

                pch = strtok(line, " ,");
                while (pch != NULL) {
                    if (str_count == obs_timestamp)
                        Z[i] = atof(pch);
                    pch = strtok(NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }

    //Copy Z to Zcpy
    MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descZ, MORSE_descZcpy);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");


    return 0;

}

int MORSE_ML_GenerateZVec_Tile_Async(char *obs_dir, int obs_timestamp) {

    //Initialization
    double flops = 0;

    //Set C diagonal equals to 1 (simple case)
    if (async == 0)

        MORSE_dlaset_Tile_Async(MorseUpperLower, 0, 1, MORSE_descC, msequence, &mrequest[1]);


    if (verbose == 1)
        fprintf(stderr, "Initializing Co-variance Matrix ...\n");


    //Generate co-variance matrix C
    MORSE_MLE_GenCovMat_Tile_Async(MORSE_descC, msequence, &mrequest[1], THETA);


    if (verbose == 1)
        fprintf(stderr, "Done ...\n");

    if (strcmp(obs_dir, "") == 0) {
        if (verbose == 1)
            fprintf(stderr,
                    "Cholesky factorization of Sigma (Generation Phase) ...\n");

        //Cholesky factorization
        int success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descC, msequence, &mrequest[1]);

        if (success != PLASMA_SUCCESS) {
            printf("Factorization cannot be performed..\n"
                   "The matrix is not positive definite\n\n");
            exit(0);
        }

        flops = flops + FLOPS_DPOTRF(N);
        if (verbose == 1)
            fprintf(stderr, " Done.\n");


        // Generate Z Vector
        MORSE_MLE_GenZVec_Tile_Async(MORSE_descZ, msequence, &mrequest[1]);

        MORSE_dtrmm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ,
                               msequence, &mrequest[1]);
    } else {
        if (verbose == 1)
            fprintf(stderr,
                    "Reading observations from dir ...\n");

        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        char *pch;
        int str_count = 0;
        char *zeroes;
        char integer_string[32];
        int i = 0;

        char *obs_dir_name = (char *) malloc(50 * sizeof(char));

        for (i = 0; i < 103; i++) {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i + 1));

            if (i < 9)
                zeroes = "00";
            else if (i >= 9 && i < 99)
                zeroes = "0";
            else
                zeroes = "";

            strcat(obs_dir_name, "/ppt.complete.Y");

            strcat(obs_dir_name, zeroes);
            strcat(obs_dir_name, integer_string);


            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count = 0;

                pch = strtok(line, " ,");
                while (pch != NULL) {
                    if (str_count == obs_timestamp)
                        Z[i] = atof(pch);
                    pch = strtok(NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }

    //Copy Z to Zcpy
    MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descZ, MORSE_descZcpy, msequence, &mrequest[1]);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");


    return 0;


}


double PLASMA_MLE_Tile(unsigned n, const double *theta, double *grad, void *my_func_data) {


//printf("******************8hi\n");
//exit(0);
    //Initialization
    double theta_hat, det = 1.0, logdet = 0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0;
    double flops = 0.0;


    PLASMA_dlaset_Tile(PlasmaUpperLower, 0, 1, descC);

    //Generate new co-variance matrix C based on new theta
    PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], theta);


    //re-store old Z
    PLASMA_dlacpy_Tile(PlasmaUpperLower, descZcpy, descZ);


    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = PLASMA_dpotrf_Tile(PlasmaLower, descC);

    STOP_TIMING(time_facto);

    if (success != PLASMA_SUCCESS) {
        printf("Factorization cannot be performed..\n"
               "The matrix is not positive definite\n\n");

        exit(0);
    }


    flops = flops + FLOPS_DPOTRF(N);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);


    START_TIMING(logdet_calculate);


    PLASMA_MLE_det_Tile_Async(descC, sequence, &request[1], &det);


    logdet = det == 0 ? 0 : log(det * det);


    STOP_TIMING(logdet_calculate);

    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        printf("log-determinant=%f\n\n", logdet);


    //Solving Linear System (L*X=Z)--->inv(L)*Z
    if (verbose == 1)
        fprintf(stderr, "Solving the linear system ...\n");


    START_TIMING(time_solve);

    //Compute triangular solve LC*X = Z
    PLASMA_dtrsm_Tile(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1, descC, descZ);

    STOP_TIMING(time_solve);

    flops = flops + FLOPS_DTRSM(PlasmaLeft, N, NRHS);

    if (verbose == 1)
        fprintf(stderr, "Calculating the log likelihood ...");
    theta_hat = -0.5 * cblas_ddot(N, Z, 1, Z, 1) - 0.5 * logdet
                - (double) (N / 2) * log(2 * PI);


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    // Print Iteration Summary
    printf("***************************************************\n");
    //printf("------logdet: %2.6f ",logdet);
//	printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
//	printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1],
           theta[2], theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Total Time: %6.2f\n", time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + logdet_calculate + time_solve));
    printf("***************************************************\n");


    iter_count++;
    return theta_hat;
}

double PLASMA_MLE_Tile_Async(unsigned n, const double *theta, double *grad, void *my_func_data) {

    //Initialization
    double theta_hat, det = 1.0, logdet = 0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0;
    double flops = 0.0;


    PLASMA_dlaset_Tile_Async(PlasmaUpperLower, 0, 1, descC, sequence, &request[1]);


    //Generate new co-variance matrix C based on new theta
    PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], theta);


    //re-store old Z
    PLASMA_dlacpy_Tile_Async(PlasmaUpperLower, descZcpy, descZ, sequence, &request[1]);



    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = PLASMA_dpotrf_Tile_Async(PlasmaLower, descC, sequence, &request[1]);


    STOP_TIMING(time_facto);

    if (success != PLASMA_SUCCESS) {
        printf("Factorization cannot be performed..\n"
               "The matrix is not positive definite\n\n");

        exit(0);
    }


    flops = flops + FLOPS_DPOTRF(N);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    //***************************************************************************************************************Stop


    //Solving Linear System (L*X=Z)--->inv(L)*Z
    if (verbose == 1)
        fprintf(stderr, "Solving the linear system ...\n");


    START_TIMING(time_solve);

    //Compute triangular solve LC*X = Z
    PLASMA_dtrsm_Tile_Async(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 1, descC, descZ, sequence,
                            &request[1]);

    STOP_TIMING(time_solve);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");

    flops = flops + FLOPS_DTRSM(PlasmaLeft, N, NRHS);

    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);


    PLASMA_MLE_det_Tile_Async(descC, sequence, &request[1], &det);


    logdet = det == 0 ? 0 : log(det * det);


    STOP_TIMING(logdet_calculate);

    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        printf("log-determinant=%f\n\n", logdet);


    if (verbose == 1)
        fprintf(stderr, "Calculating the log likelihood ...");
    theta_hat = -0.5 * cblas_ddot(N, Z, 1, Z, 1) - 0.5 * logdet
                - (double) (N / 2) * log(2 * PI);


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    // Print Iteration Summary
    printf("***************************************************\n");
    //printf("------logdet: %2.6f ",logdet);
    //printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    //printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1],
           theta[2], theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Total Time: %6.2f\n", time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + logdet_calculate + time_solve));
    printf("***************************************************\n");


    iter_count++;
    return theta_hat;
}


double MORSE_MLE_Tile(unsigned n, const double *theta, double *grad, void *my_func_data) {

    //Initialization
    double theta_hat, det = 1.0, logdet = 0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0;
    double flops = 0.0;


    MORSE_dlaset_Tile(MorseUpperLower, 0, 1, MORSE_descC);

    //Generate new co-variance matrix C based on new theta
    MORSE_MLE_GenCovMat_Tile_Async(MORSE_descC, msequence, &mrequest[1], theta);


    //re-store old Z
    MORSE_dlacpy_Tile(MorseUpperLower, MORSE_descZcpy, MORSE_descZ);


    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = MORSE_dpotrf_Tile(MorseLower, MORSE_descC);

    STOP_TIMING(time_facto);

    if (success != PLASMA_SUCCESS) {
        printf("Factorization cannot be performed..\n"
               "The matrix is not positive definite\n\n");

        exit(0);
    }


    flops = flops + FLOPS_DPOTRF(N);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);


    START_TIMING(logdet_calculate);


    MORSE_MLE_det_Tile_Async(MORSE_descC, msequence, &mrequest[1], &det);


    logdet = det == 0 ? 0 : 2 * log(det);


    STOP_TIMING(logdet_calculate);

    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        printf("log-determinant=%f\n\n", logdet);


    //Solving Linear System (L*X=Z)--->inv(L)*Z
    if (verbose == 1)
        fprintf(stderr, "Solving the linear system ...\n");


    START_TIMING(time_solve);

    //Compute triangular solve LC*X = Z
    MORSE_dtrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ);

    STOP_TIMING(time_solve);

    flops = flops + FLOPS_DTRSM(MorseLeft, N, NRHS);

    if (verbose == 1)
        fprintf(stderr, "Calculating the log likelihood ...");
    theta_hat = -0.5 * cblas_ddot(N, Z, 1, Z, 1) - 0.5 * logdet
                - (double) (N / 2) * log(2 * PI);


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    // Print Iteration Summary
    printf("***************************************************\n");
    //printf("------logdet: %2.6f ",logdet);
    //printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    //printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1],
           theta[2], theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Total Time: %6.2f\n", time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + logdet_calculate + time_solve));
    printf("***************************************************\n");


    iter_count++;
    return theta_hat;
}

double MORSE_MLE_Tile_Async(unsigned n, const double *theta, double *grad, void *my_func_data) {

    //Initialization
    double theta_hat, det = 1.0, logdet = 0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0;
    double flops = 0.0;


    MORSE_dlaset_Tile_Async(MorseUpperLower, 0, 1, MORSE_descC, msequence, &mrequest[1]);


    //Generate new co-variance matrix C based on new theta
    MORSE_MLE_GenCovMat_Tile_Async(MORSE_descC, msequence, &mrequest[1], theta);


    //re-store old Z
    MORSE_dlacpy_Tile_Async(MorseUpperLower, MORSE_descZcpy, MORSE_descZ, msequence, &mrequest[1]);



    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = MORSE_dpotrf_Tile_Async(MorseLower, MORSE_descC, msequence, &mrequest[1]);


    STOP_TIMING(time_facto);

    if (success != MORSE_SUCCESS) {
        printf("Factorization cannot be performed..\n"
               "The matrix is not positive definite\n\n");

        exit(0);
    }


    flops = flops + FLOPS_DPOTRF(N);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    //***************************************************************************************************************Stop


    //Solving Linear System (L*X=Z)--->inv(L)*Z
    if (verbose == 1)
        fprintf(stderr, "Solving the linear system ...\n");


    START_TIMING(time_solve);

    //Compute triangular solve LC*X = Z
    MORSE_dtrsm_Tile_Async(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, 1, MORSE_descC, MORSE_descZ, msequence,
                           &mrequest[1]);

    STOP_TIMING(time_solve);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");

    flops = flops + FLOPS_DTRSM(MorseLeft, N, NRHS);

    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);


    MORSE_MLE_det_Tile_Async(MORSE_descC, msequence, &mrequest[1], &det);


    logdet = det == 0 ? 0 : log(det * det);


    STOP_TIMING(logdet_calculate);

    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    if (verbose == 1)
        printf("log-determinant=%f\n\n", logdet);


    if (verbose == 1)
        fprintf(stderr, "Calculating the log likelihood ...");
    theta_hat = -0.5 * cblas_ddot(N, Z, 1, Z, 1) - 0.5 * logdet
                - (double) (N / 2) * log(2 * PI);


    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    // Print Iteration Summary
    printf("***************************************************\n");
    //printf("------logdet: %2.6f ",logdet);
    //printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    //printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1],
           theta[2], theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Total Time: %6.2f\n", time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + logdet_calculate + time_solve));
    printf("***************************************************\n");


    iter_count++;
    return theta_hat;
}

void ML_finalize() {

    //writetofile("theta_hat.txt", iter, 1, theta_hat_plot,		sizeof(theta) / sizeof(double));
    //writetofile("theta.txt", iter, 1, theta, sizeof(theta) / sizeof(double));
    free(C);
    free(Z);
    free(Zcpy);
    if (check == 1) {
        free(Ccpy);

    }
}


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


int GenerateXYLoc(int n, char *locs_file) {


    int i;

    if (strcmp(locs_file, "") == 0) {
        //Uniform random generation of distance matrix (For testing Phase)
        if (verbose == 1)
            fprintf(stderr, "Initializing Locations ...\n");

        for (i = 0; i < n; i++) {
            X[i] = (R - 0.5 + uniform_distribution(MIN_RAND, MAX_RAND))
                   * sqrt((double) n);
            Y[i] = (L - 0.5 + uniform_distribution(MIN_RAND, MAX_RAND))
                   * sqrt((double) n);
        }
    } else {
        if (verbose == 1)
            fprintf(stderr, "Reading Locations from file ...\n");


        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        char *pch;
        int str_count = 0;

        fp = fopen(locs_file, "r");
        if (fp == NULL)
            exit(EXIT_FAILURE);

        //avoid header
        read = getline(&line, &len, fp);
        while ((read = getline(&line, &len, fp)) != -1) {
            str_count = 0;

            pch = strtok(line, " ,");
            while (pch != NULL) {
                if (str_count == 1)
                    X[i] = atof(pch);
                else if (str_count == 2)
                    Y[i] = atof(pch);
                pch = strtok(NULL, " ,");
                str_count++;
            }

            //	printf("X[i]=: %f, Y[i]=: %f\n", X[i],Y[i]);
            i++;

        }

        fclose(fp);
        if (line)
            free(line);
    }

    return 0;


}


int PLASMA_MLE_GenCovMat_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
                                    PLASMA_request *request, double *theta) {

    //Dynamic scheduler functions
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m, n;
    int ldam;
    int tempnn, tempmm;
    PLASMA_desc A = *descA;


    //printf("%f, %f, %f, \n", theta[0],theta[1],theta[2]);
    //exit(0);

    //get Quark for PLASMA
    PLASMA_Get_Quark(&quark);

    //Set various task level flags. This flag data structure is then provided when the task is created/inserted. Each flag can take a value which is either an integer or a pointer.
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE,
                        (intptr_t) sequence->quark_sequence);

    //mt is the number of tile rows of the sub-matrix -- nt is the number of tile columns of the sub-matrix -- mb the number of rows in a tile
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A.nt; n++) {

            tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
            QUARK_CORE_dGenCovMat(quark, &task_flags, tempmm, tempnn, A(m, n), ldam, m * A.mb, n * A.nb, theta);
        }

    }

    return PLASMA_SUCCESS;
}


int PLASMA_MLE_GenZVec_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
                                  PLASMA_request *request) {


    //Dynamic scheduler functions
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m;
    int ldam;
    int tempmm;
    PLASMA_desc A = *descA;

    //get Quark for PLASMA
    PLASMA_Get_Quark(&quark);

    //Set various task level flags. This flag data structure is then provided when the task is created/inserted. Each flag can take a value which is either an integer or a pointer.
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE,
                        (intptr_t) sequence->quark_sequence);


    double *r = (double *) malloc(A.m * sizeof(double));

    //Uniform Random Generator
    int *iseed = (int *) malloc(4 * sizeof(int));
    iseed[0] = 371;
    iseed[1] = 371;
    iseed[2] = 371;
    iseed[3] = 371;

    //Uniform random generation of e --  ei~N(0,1)
    LAPACKE_dlarnv(3, iseed, A.m, r);


    //mt is the number of tile rows of the sub-matrix -- nt is the number of tile columns of the sub-matrix -- mb the number of rows in a tile
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = BLKLDD(A, m);


        QUARK_CORE_dGenZVec(quark, &task_flags, tempmm, 1, A(m, 0), ldam, m * A.mb, 0, r);
    }

    return PLASMA_SUCCESS;
}


int PLASMA_MLE_det_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
                              PLASMA_request *request, double *det_acc) {


    //Dynamic scheduler functions
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m, n;
    int ldam;
    int tempnn, tempmm;
    PLASMA_desc A = *descA;

    *det_acc = 1;


    //get Quark for PLASMA
    PLASMA_Get_Quark(&quark);

    //Set various task level flags. This flag data structure is then provided when the task is created/inserted. Each flag can take a value which is either an integer or a pointer.
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE,
                        (intptr_t) sequence->quark_sequence);


    //mt is the number of tile rows of the sub-matrix -- nt is the number of tile columns of the sub-matrix -- mb the number of rows in a tile
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A.nt; n++) {
            //generate the Lower and diagonal tiles if symmetric
            // if(part == 1 && n > m) break;
            tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
            QUARK_CORE_det(quark, &task_flags, tempmm, tempnn, A(m, n), ldam, m * A.mb, n * A.nb, det_acc);


        }

    }


    QUARK_Barrier(quark);

    return PLASMA_SUCCESS;
}


void QUARK_CORE_dGenCovMat(Quark *quark, Quark_Task_Flags *task_flags, int m,
                           int n, double *A, int lda, int m0, int n0, double *theta) {
    /*DAG_CORE_GenCovMat;
     Called by the master thread. Add a new task to the scheduler, providing the data pointers, sizes, and dependency information.
     This function provides the main user interface for the user to write data-dependent algorithms.*/
//printf("%f, %f, %f, \n", theta[0],theta[1],theta[2]);
//exit(0);


    QUARK_Insert_Task(quark, CORE_dGenCovMat_quark, task_flags, sizeof(int), &m,
                      VALUE, sizeof(int), &n, VALUE, sizeof(double) * lda * n, A,
                      OUTPUT | LOCALITY, sizeof(int), &lda, VALUE, sizeof(int), &m0,
                      VALUE, sizeof(int), &n0, VALUE, 3 * sizeof(double), theta, INPUT, 0);

}

void QUARK_CORE_dGenZVec(Quark *quark, Quark_Task_Flags *task_flags, int m,
                         int n, double *A, int lda, int m0, int n0, double *r) {
    /*DAG_CORE_GenCovMat;
     Called by the master thread. Add a new task to the scheduler, providing the data pointers, sizes, and dependency information.
     This function provides the main user interface for the user to write data-dependent algorithms.*/
    QUARK_Insert_Task(quark, CORE_dGenZVec_quark, task_flags, sizeof(int), &m,
                      VALUE, sizeof(int), &n, VALUE, sizeof(double) * lda * n, A,
                      OUTPUT | LOCALITY, sizeof(int), &lda, VALUE, sizeof(int), &m0,
                      VALUE, sizeof(int), &n0, VALUE, sizeof(double) * lda * n, r, INPUT, 0);


}


void QUARK_CORE_det(Quark *quark, Quark_Task_Flags *task_flags, int m,
                    int n, double *A, int lda, int m0, int n0, double *determinant) {
    /*DAG_CORE_GenCovMat;
     Called by the master thread. Add a new task to the scheduler, providing the data pointers, sizes, and dependency information.
     This function provides the main user interface for the user to write data-dependent algorithms.*/
    QUARK_Insert_Task(quark, CORE_det_quark, task_flags, sizeof(int), &m,
                      VALUE, sizeof(int), &n, VALUE, sizeof(double) * lda * n, A,
                      INPUT, sizeof(int), &lda, VALUE, sizeof(int), &m0,
                      VALUE, sizeof(int), &n0, VALUE, sizeof(double), determinant, INOUT, 0);

}

void CORE_dGenCovMat_quark(Quark *quark) {

    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double *theta;


    quark_unpack_args_7(quark, m, n, A, lda, m0, n0, theta);

//printf("%f, %f, %f, \n", theta[0],theta[1],theta[2]);
//exit(0);

    matcov_comp_Tile(A, m, n, m0, n0, theta);

}

void CORE_dGenZVec_quark(Quark *quark) {

    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double *r;

    quark_unpack_args_7(quark, m, n, A, lda, m0, n0, r);
    veccov_comp_Tile(A, m, n, m0, n0, r);

}


void CORE_det_quark(Quark *quark) {

    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double *determinant;

    quark_unpack_args_7(quark, m, n, A, lda, m0, n0, determinant);
    det_comp_Tile(A, m, n, m0, n0, determinant);
}

void matcov_comp_Tile(double *A, int m, int n, int m0, int n0, double *theta) {

    int i, j;
    int x = m0;
    int y = n0;

    double dist = 0.0;
    double expr = 0.0;
    double con = 0.0;
//printf("%f\n", tgamma(2.5));
//exit(0);
    //printf("result: %f\n",bessk(3,5));//jn(n,x)
    //exit(0);

    if (strcmp(kernel, "MATERN_KERENL") == 0) {
        con = pow(2, (theta[2] - 1)) * tgamma(theta[2]);
        con = 1.0 / con;
        con = theta[0] * con;

        for (i = 0; i < m; i++) {
            y = n0;
            for (j = 0; j < n; j++) {


//			if(A[i+j*m]!=1)
//			{
                dist = CalculateDistance(X[x], Y[x], X[y], Y[y]);
                expr = dist / theta[1];
                if (expr == 0)
                    expr = 1e-10;


                A[i + j * m] = con * pow(expr, theta[2]) * bessk(theta[2], expr); // Matern Function
//			}



                y++;
            }
            x++;
        }
    } else {
        printf("start from: %d , to: %d \n", m, ((m - 1) + n * m));
        for (i = m0; i < m; i++) {
            y = n0;
            for (j = n0; j < n; j++) {
                A[i + j * m] =
                        A[i + j * m] == 1 ?
                        1 :
                        exp(-CalculateDistance(X[x], Y[x], X[y], Y[y]) / theta[0]);

                y++;
            }
            x++;
        }

    }
}


void det_comp_Tile(double *A, int m, int n, int m0, int n0, double *res) {

//*res=1;
//printf("yes\n");
    int i, j;
    int x = m0;
    int y = n0;

    for (i = 0; i < m; i++) {
        y = n0;
        for (j = 0; j < n; j++) {
            if (x == y) {
                *res *= A[i + j * m];
//printf("A[i+j*m: %f \n",A[i+j*m]);
            }

            y++;
        }
        x++;
    }

//printf("*res= %f\n", *res);
}


void veccov_comp_Tile(double *Z, int m, int n, int m0, int n0, double *r) {

    int i = 0;

    for (i = 0; i < m; i++)
        Z[i] = r[m0 + i];
}

void print_matrix(char *desc, int m, int n, double *a, int lda) {
    int i, j;
    printf("\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf(" %6.4f", a[i + lda * j]);
        printf("\n");
    }
}

//*******************************************************************************(1)
void CORE_dGenCovMat_starpu(void *buffers[], void *cl_arg) {
    int m, n, lda, m0, n0;
    double *theta;
    double *A;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &A, &lda, &m0, &n0, &theta);
    matcov_comp_Tile(A, m, n, m0, n0, theta);
}

static struct starpu_codelet cl_dGenCovMat =
        {
                .where = STARPU_CPU,
                .cpu_funcs = {CORE_dGenCovMat_starpu},
                .nbuffers = 0
        };


int MORSE_MLE_GenCovMat_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request,
                                   double *theta) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m, n, m0, n0;
    int ldam;
    int tempmm, tempnn;
    MORSE_desc_t A = *descA;

    struct starpu_codelet *cl = &cl_dGenCovMat;


    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = mBLKLDD(descA, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;

            double *data = (double *) morse_getaddr_ccrb(descA, m, n);

            m0 = m * A.mb;
            n0 = n * A.nb;
            //      printf("m0=%d, ",m0);
            //	printf("n0=%d, \n",n0);

            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &data, sizeof(double *),
                               STARPU_VALUE, &ldam, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_VALUE, &theta, sizeof(double *),
                               0);
        }
    }

    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
    return MORSE_SUCCESS;
}

//****************************************************************************(2)
void CORE_GenZVec_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double *r;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &A, &lda, &m0, &n0, &r);
    veccov_comp_Tile(A, m, n, m0, n0, r);
}

static struct starpu_codelet cl_GenZVec =
        {
                .where = STARPU_CPU,
                .cpu_funcs = {CORE_GenZVec_starpu},
                .nbuffers = 0
        };


int MORSE_MLE_GenZVec_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m, m0, n0;
    int ldam;
    int tempmm, tempnn;
    MORSE_desc_t A = *descA;


    double *r = (double *) malloc(A.m * sizeof(double));

    //Uniform Random Generator
    int *iseed = (int *) malloc(4 * sizeof(int));
    iseed[0] = 371;
    iseed[1] = 371;
    iseed[2] = 371;
    iseed[3] = 371;

    //Uniform ranodm generation of e -- ei~N(0, 1)
    LAPACKE_dlarnv(3, iseed, A.m, r);

    struct starpu_codelet *cl = &cl_GenZVec;


    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = mBLKLDD(descA, m);

        tempnn = 1;

        double *data = (double *) morse_getaddr_ccrb(descA, m, 0);

        m0 = m * A.mb;
        n0 = 0;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempnn, sizeof(int),
                           STARPU_VALUE, &data, sizeof(double *),
                           STARPU_VALUE, &ldam, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_VALUE, &r, sizeof(double),
                           0);

    }

    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
    return MORSE_SUCCESS;
}


//*****************************************************************************(3)
void CORE_det_starpu(void *buffers[], void *cl_arg) {
    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double *determinant;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &A, &lda, &m0, &n0, &determinant);
    det_comp_Tile(A, m, n, m0, n0, determinant);

//printf("det=%f:\n", determinant);
}

static struct starpu_codelet cl_det =
        {
                .where = STARPU_CPU,
                .cpu_funcs = {CORE_det_starpu},
                .nbuffers = 0
        };


int
MORSE_MLE_det_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request, double *det_acc) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, morse, sequence, request);

    int m, n, m0, n0;
    int ldam;
    int tempmm, tempnn;
    MORSE_desc_t A = *descA;
    *det_acc = 1.0;
//double det_val=*det_acc;
    struct starpu_codelet *cl = &cl_det;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = mBLKLDD(descA, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;

            double *data = (double *) morse_getaddr_ccrb(descA, m, n);

            m0 = m * A.mb;
            n0 = n * A.nb;

            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &data, sizeof(double *),
                               STARPU_VALUE, &ldam, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_VALUE, &det_acc, sizeof(double *),
                               0);
//printf("det_val: %f\n",det_val);
        }
    }
    RUNTIME_barrier(morse);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
    return MORSE_SUCCESS;
}

int countlines(char *filename) {
    // count the number of lines in the file called filename
    FILE *fp = fopen(filename, "r");
    int ch = 0;
    int lines = 0;

    if (fp == NULL) {
        printf("cannot open locations file\n");
        return 0;
    }


    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == '\n')
            lines++;
    }


    fclose(fp);

    //Excluding header line
    return (lines - 1);
}


static double bessi0(double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
    double ax, ans;
    double y;


    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75, y = y * y;
        ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                                           + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    } else {
        y = 3.75 / ax;
        ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                                                        + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                                                                                                      + y *
                                                                                                        (-0.2057706e-1 +
                                                                                                         y *
                                                                                                         (0.2635537e-1 +
                                                                                                          y *
                                                                                                          (-0.1647633e-1
                                                                                                           + y *
                                                                                                             0.392377e-2))))))));
    }
    return ans;
}


static double bessi1(double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
    double ax, ans;
    double y;


    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75, y = y * y;
        ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                                                                   + y * (0.2658733e-1 +
                                                                          y * (0.301532e-2 + y * 0.32411e-3))))));
    } else {
        y = 3.75 / ax;
        ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
                                                       - y * 0.420059e-2));
        ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
                                                     + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

static double bessk0(double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0.  */
/*------------------------------------------------------------*/
{
    double y, ans;

    if (x <= 2.0) {
        y = x * x / 4.0;
        ans = (-log(x / 2.0) * bessi0(x)) + (-0.57721566 + y * (0.42278420
                                                                + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
                                                                                                             + y *
                                                                                                               (0.10750e-3 +
                                                                                                                y *
                                                                                                                0.74e-5))))));
    } else {
        y = 2.0 / x;
        ans = (exp(-x) / sqrt(x)) * (1.25331414 + y * (-0.7832358e-1
                                                       + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
                                                                                                       + y *
                                                                                                         (-0.251540e-2 +
                                                                                                          y *
                                                                                                          0.53208e-3))))));
    }
    return ans;
}


static double bessk1(double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
    double y, ans;

    if (x <= 2.0) {
        y = x * x / 4.0;
        ans = (log(x / 2.0) * bessi1(x)) + (1.0 / x) * (1.0 + y * (0.15443144
                                                                   + y * (-0.67278579 +
                                                                          y * (-0.18156897 + y * (-0.1919402e-1
                                                                                                  + y * (-0.110404e-2 +
                                                                                                         y *
                                                                                                         (-0.4686e-4)))))));
    } else {
        y = 2.0 / x;
        ans = (exp(-x) / sqrt(x)) * (1.25331414 + y * (0.23498619
                                                       + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
                                                                                                       + y *
                                                                                                         (0.325614e-2 +
                                                                                                          y *
                                                                                                          (-0.68245e-3)))))));
    }
    return ans;
}


double bessk(int n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n >= 0*/
/* Note that for x == 0 the functions bessy and bessk are not */
/* defined and a blank is returned.                           */
/*------------------------------------------------------------*/
{
    int j;
    double bk, bkm, bkp, tox;


    if (n < 0 || x == 0.0) {
        double dblank;
        //   setdblank_c( &dblank );
        return (dblank);
    }
    if (n == 0)
        return (bessk0(x));
    if (n == 1)
        return (bessk1(x));

    tox = 2.0 / x;
    bkm = bessk0(x);
    bk = bessk1(x);
    for (j = 1; j < n; j++) {
        bkp = bkm + j * tox * bk;
        bkm = bk;
        bk = bkp;
    }
    return bk;
}
