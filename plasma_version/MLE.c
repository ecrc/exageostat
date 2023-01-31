/**
 *
 * @file MLE.c
 *  
 *  MLE is a software package provided by KAUST,
 *  King Abdullah University of Science and Technology - ECRC
 *
 * @version 1.0.0
 * @author Sameh Abdulsh
 * @date 2021-01-24
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <nlopt.h>
#include <math.h>
#include <plasma.h>
#include <chameleon.h>
#include <starpu.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "flops.h"
#include <argp.h>
#include <gsl/gsl_sf_bessel.h>


// For Bessel Function (Do not change)
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
//tiling
#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )
#define mBLKLDD(A, k) A->get_blkldd( A,k )
#define A(m,n) (double *)plasma_getaddr(A, m, n)
//Sameh Testing MOde
#define MIN_RAND       -0.4
#define MAX_RAND        0.4
#define R        2
#define L       3
//#define THETA1 1 //weak=0.03, medium=0.1,  strong=0.3
//#define THETA2 0.05
//#define THETA3 1
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
static void print_matrix(char* desc, int m, int n, double* a, int lda);
int PLASMA_MLE_GenZVec_Tile_Async(PLASMA_desc *descA, int * r, PLASMA_sequence *sequence,
        PLASMA_request *request);
int EXAGEOSTAT_MLE_GenZVec_Tile_Async(CHAM_desc_t *descA, int * r,RUNTIME_sequence_t *sequence,
        RUNTIME_request_t *request);
static int GenerateXYLoc(int n, char * locs_file,int based_sys);
static double CalculateDistance(double x1, double y1, double x2, double y2);
static double uniform_distribution(double rangeLow, double rangeHigh);
//tiling
static void CORE_dGenCovMat_quark(Quark *quark);
static void CORE_dGenZVec_quark(Quark *quark);
static int PLASMA_MLE_GenCovMat_Tile_Async(PLASMA_desc *descA,
        PLASMA_sequence *sequence, PLASMA_request *request, double * theta);
static int EXAGEOSTAT_MLE_GenCovMat_Tile_Async(CHAM_desc_t *descA,
        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request, double * theta);
static void QUARK_CORE_dGenCovMat(Quark *quark, Quark_Task_Flags *task_flags,
        int m, int n, double *A, int lda, int m0, int n0, double * theta);
static void matcov_comp_Tile(double * A, int m, int n, int m0, int n0, double * theta);
static void  CORE_det_quark(Quark *quark);
static void  det_comp_Tile(double * A, int m, int n, int m0, int n0,double * determinant) ;
static void QUARK_CORE_det(Quark *quark, Quark_Task_Flags *task_flags, int m,
        int n, double *A, int lda, int m0, int n0, double * determinant );
static int PLASMA_MLE_det_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
        PLASMA_request *request,double * result);
static int PLASMA_MLE_det_Tile(PLASMA_desc *descA, PLASMA_sequence *sequence,
        PLASMA_request *request,double * det_acc);
static int EXAGEOSTAT_MLE_det_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence,
        RUNTIME_request_t *request,double * result);
void QUARK_CORE_dGenZVec(Quark *quark, Quark_Task_Flags *task_flags, int m,
        int n, double *A, int lda, int m0, int n0, double * r);
void veccov_comp_Tile(double * Z, int m, int n, int m0, int n0, double *r);
int PLASMA_ML_GenerateZVec_Tile(char * obs_di,int obs_timestamp, int * r);
double PLASMA_MLE_Tile(unsigned n, const double * theta, double * grad, void * my_func_data);
void ML_finalize();
int PLASMA_ML_GenerateZVec_Tile_Async(char * obs_dir, int obs_timestamp,int * r);
double PLASMA_MLE_Tile_Async(unsigned n, const double * theta, double * grad, void * my_func_data);
int EXAGEOSTAT_ML_GenerateZVec_Tile(char * obs_dir, int obs_timestamp, int * r);
double EXAGEOSTAT_MLE_Tile(unsigned n, const double * theta, double * grad, void * my_func_data);
int EXAGEOSTAT_ML_GenerateZVec_Tile_Async(char * obs_dir, int obs_timestamp, int * seeds);
double EXAGEOSTAT_MLE_Tile_Async(unsigned n, const double * theta, double * grad, void * my_func_data);
int EXAGEOSTAT_MLE_det_Tile(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence,
        RUNTIME_request_t *request,double * det_acc) ;
int countlines(char *filename);
double bessy( int n, double x );
static double bessy1( double x );
static double bessy0( double x );
static double bessj1( double x );
static double bessj0( double x );
void theta_parser(double * theta_vec,char * kern);
//for experiments
void write_to_file(char * path, int matrix_size,int ncores,int tile_size, int test, char * ikernel, int based_system, int async, char *obs_dir,int obs_timestamp,double total_exec_time,double avg_exec_time_per_iter, double avg_flops_per_iter);



double * C ;
double * Z;
double * Zcpy;
double * Ccpy;
//********************************************************************
PLASMA_desc *descC = NULL;
PLASMA_desc *descZ = NULL;
PLASMA_desc *descZcpy = NULL;
CHAM_desc_t *CHAM_descC = NULL;
CHAM_desc_t *CHAMELEON_descZ = NULL;
CHAM_desc_t *CHAMELEON_descZcpy = NULL;
double * X;
double * Y;
//PLASMA sequence uniquely identifies a set of asynchronous function calls sharing common exception handling.
PLASMA_sequence *sequence;
//PLASMA request uniquely identifies each asynchronous function call.
PLASMA_request request[19] = { PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER,
    PLASMA_REQUEST_INITIALIZER, PLASMA_REQUEST_INITIALIZER };
//MORSE sequence uniquely identifies a set of asynchronous function calls sharing common exception handling.
RUNTIME_sequence_t *msequence;
//MORSE request uniquely identifies each asynchronous function call.
RUNTIME_request_t mrequest[19] = { CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS,
    CHAMELEON_SUCCESS, CHAMELEON_SUCCESS };
//***************************************************************************************

int N,NRHS,LDC,LDZ,ts,async,check,verbose;
int iter_count=0;
int based_sys;
double *target_theta;
double *initial_theta;
char *kernel;
char *ikernel;

//to collect results
double avg_exec_time_per_iter;
double total_exec_time;
double avg_flops_per_iter;
double final_theta_hat;

const char *argp_program_version ="Version 1.0";
const char *argp_program_bug_address ="<sameh.abdulah@kaust.edu.sa>";
struct arguments
{
    //  char *args[0];          /* ncores, N, ts, locs_file, obs_dir, and obs_timestamp */
    int test;                   /* The -t flag (USE Synthetic Dataset*/
    int check ;                /* The -c flag */
    char *zvecs ;                /* The number of Z vectors to be tested  */
    int verbose;             /* The -v flag */
    char *ncores;              //Number of cores
    char *N;                   //Matrix size (Only in the case of Synthetic dataset)
    char *ts;                  //Tiling size
    char *kernel;
    char *ikernel;
    int chameleon;		  // 0--> plasma 1-->chameleon  The -sys flag
    int async; 		  // 0-->tile  1-->tile_async       The -async flag
    char *locs_file;        // Locations files (in the case of real dataset
    char *obs_dir;          // Observations directory in the case of real dataset
    int timestamp;       // Observations timestamp

};

static struct argp_option options[] =
{
    {"test", 't', 0, 0, "Execute in test mode"},
    {"check", 'c', 0, 0, "Produce check output"},
    {"zvecs",   'z', "ZVECS", 0, "number of Z vectors to be tested"},
    {"verbose", 'v', 0, 0, "Produce verbose output"},
    {"ncores",   'n', "NCORES", 0, "Number of Cores"},
    {"N",   's', "MATRIX_SIZE", 0, "Synthetic Matrix Size"},
    {"ts",   'e', "TILE_SIZE", 0, "Number of Tiles"},
    {"kernel",   'k', "KERNEL", 0, "Computation Model"},
    {"ikernel",   'i', "IKERNEL", 0, "Initial theta(s) used for Testing case Generation"},
    {"chameleon",   'b', 0, 0, "Use Chameleon Instead of Plasma"},
    {"async", 'a', 0, 0, "Asynchronous"},
    {"locs_file",   'l', "LOCATIONS_FILE", 0, "Read Locations from this Location File"},
    {"obs_dir",  'o', "OBSERVATIONS_DIRECTORY", 0, "Read Observations from this Directory Path"},
    {"timestamp",  'p', "TIMESTAMP", 0, "Observation Timestamp"},
    {0}
};



    static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
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
        case 'i':
            arguments->ikernel = arg;
            break;
        case 'z':
            arguments->zvecs = arg;
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



int main(int argc, char **argv) {

    //initialization
    double *starting_theta;
    char *token;
    int i = 0;
    int locs;
    double time_opt = 0.0;
    double max_theta_hat = 0.0;
    int zvecs = 1;
    int iseed = 0;
    //For experiments
    avg_exec_time_per_iter = 0;
    total_exec_time = 0;
    avg_flops_per_iter = 0;


    struct arguments arguments;
    FILE *outstream;
    arguments.test = 0;                  /* execute in test mode */
    arguments.check = 0;                 /* The -c flag */
    arguments.verbose = 0;             /* The -v flag */
    arguments.zvecs = "1";             /* The number of Z vectors to be tested*/
    arguments.chameleon = 0;           // 0--> plasma 1-->chameleon  The -sys flag      
    arguments.async = 0;               // 0-->tile  1-->tile_async       The -async flag
    arguments.kernel = ""; 	//   target Kernel
    arguments.ikernel = ""; 	//   initial Kernel
    arguments.ncores = "1";              //Number of cores
    arguments.N = "0";                   //Matrix size (Only in the case oif sysnthetic dataset)
    arguments.ts = "0";                  //Tiling size   
    arguments.locs_file = "";        // Locations files (in the case of real dataset   
    arguments.obs_dir = "";          // Observations directory in the case of real dataset
    arguments.timestamp = 0;       // obstervations timestamp  




    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    int test = arguments.test;
    int ncores = atoi(arguments.ncores);
    N = atoi( arguments.N);
    zvecs = atoi(arguments.zvecs);
    kernel=arguments.kernel;
    ikernel=arguments.ikernel;
    ts = atoi(arguments.ts);
    based_sys = arguments.chameleon;// 0--> plasma 1-->chameleon
    async= arguments.async; // 0-->tile  1-->tile_async
    char *locs_file= arguments.locs_file;
    char *obs_dir= arguments.obs_dir;
    int obs_timestamp=arguments.timestamp;
    check = arguments.check;
    verbose = arguments.verbose;
    int num_params = 6;
    //stop gsl error handler
    gsl_set_error_handler_off () ;

    starting_theta=(double *) malloc(num_params * sizeof(double));
    double * avg_target_theta=(double *) malloc(num_params * sizeof(double));

    for (i = 0; i < num_params; i++)
    {
        avg_target_theta[i] = 6;
    }

    starting_theta[0] = 0.3;
    starting_theta[1] = 0.3;
    starting_theta[2] = 0.3;
    starting_theta[3] = 0.3;
    starting_theta[4] = 0.3;
    starting_theta[5] = 0.3;

    target_theta  = (double *) malloc(num_params * sizeof(double));
    initial_theta = (double *) malloc(num_params * sizeof(double));

    if(test ==0 && (strcmp(locs_file,"")==0 || 
                strcmp(obs_dir,"") ==0 || obs_timestamp==0 || N >0 ))
    {
        printf("\nIn real mode: please use locs_file, obs_dir, and obs_timestamp arguments only and ignore N (Matrix Size)\n\n");
        exit(0);
    }
    else if(test ==1 && (strcmp(locs_file,"")!=0 || strcmp(obs_dir,"") !=0 || obs_timestamp!=0 || N ==0 ))

    {
        printf("\nIn test mode: please use N (Matrix Size), ikernel and ignore locs_file, obs_dir, and obs_timestamp arguments\n\n");
        exit(0);

    }

    //kernel parsing
    if (test==1)
        theta_parser(initial_theta,ikernel);

    theta_parser(target_theta,kernel);


    if(strcmp(locs_file,"")!=0)
    {
        printf("locs file: %s\n",locs_file);
        locs= countlines(locs_file);
        N=locs;
    }


    NRHS = 1;
    LDC = N;
    LDZ = N;

    printf("N: %d\n",N);
    //Memory Allocation
    X = (double *) calloc((size_t)N, sizeof(double));
    Y = (double *) calloc(N , sizeof(double));
    C = (double *) calloc((size_t)LDC * (size_t)N , sizeof(double));
    Z = (double *) calloc(LDZ * NRHS , sizeof(double));
    Zcpy = (double *) calloc(LDZ * NRHS , sizeof(double));

    //printf("RLIMIT_DATA%d",RLIMIT_DATA);

    /* Check if unable to allocate memory */
    if ((!C) || (!Z) || (!Zcpy) ||(!X) || (!Y)) {
        printf("Out of Memory for C,  Zcpy , Z, X, and Y\n ");
        return -2;
    }


    //	GenerateXYLoc(N,locs_file,based_sys);

    //exit(0);
    nlopt_opt opt;
    //Initial nlopt (Optimization) 
    opt=nlopt_create(/*NLOPT_GN_CRS2_LM*/NLOPT_LN_BOBYQA , num_params);

    double lb[3]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    double up[3]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    //if(test==1)
    for(i=0;i<num_params;i++)
    {
        if(target_theta[i]!=-1)
        {
            lb[i]=target_theta[i];
            up[i]=target_theta[i];
            starting_theta[i]=target_theta[i];
        }
    }

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, up);
    double * opt_f;
    nlopt_set_xtol_rel(opt, 1e-5);

    //Uniform Random Generator
    int iseed[4]={seed, seed, seed, 1};

    double *r= (double *) malloc (N * zvecs * sizeof(double));

    //Uniform ranodm generation of e -- ei~N(0, 1)
    LAPACKE_dlarnv(3, iseed, N * zvecs, r);

    GenerateXYLoc(N, locs_file, based_sys);

    for(i=0;i<zvecs;i++)
    {
        if(based_sys==0)
        {
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

            PLASMA_Desc_Create(&descZ, Z, PlasmaRealDouble, ts, ts,ts*ts, LDZ, NRHS , 0, 0, LDZ, NRHS);

            PLASMA_Desc_Create(&descZcpy, Zcpy, PlasmaRealDouble, ts, ts, ts*ts, LDZ, NRHS, 0, 0,LDZ, NRHS);

            // Generate Observations Vector (Z) for testing phase
            if(async==0)
            {
                PLASMA_ML_GenerateZVec_Tile(obs_dir,obs_timestamp,&r[i*N]);
                START_TIMING(time_opt);

                //Maximum Likelihood (Optimization Phase using nlopt library)
                nlopt_set_max_objective(opt, PLASMA_MLE_Tile, NULL);
                nlopt_optimize(opt, starting_theta, &opt_f);

                STOP_TIMING(time_opt);
            }

            else
            {
                PLASMA_ML_GenerateZVec_Tile_Async(obs_dir,obs_timestamp,&r[i*N]);

                START_TIMING(time_opt);

                nlopt_set_max_objective(opt, PLASMA_MLE_Tile_Async, NULL);
                nlopt_optimize(opt, starting_theta, &opt_f);

                STOP_TIMING(time_opt);
            }

        }

        // MORSE based- system
        else
        {

            /*-------------------------------------------------------------
             *              Initialization of Morse
             */
            CHAMELEON_Init(ncores,1);  //CHAMELEON_init(NCPU,NGPU);
            CHAMELEON_Disable(PLASMA_AUTOTUNING);
            CHAMELEON_Set(PLASMA_TILE_SIZE, ts);
            CHAMELEON_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING);
            CHAMELEON_Enable(PLASMA_WARNINGS);
            CHAMELEON_Enable(PLASMA_ERRORS);
            /*-------------------------------------------------------------*/

            //			GenerateXYLoc(N,locs_file,based_sys);

            //Identifies a set of routines sharing common exception handling.
            CHAMELEON_Sequence_Create(&msequence);

            //Create matrix descriptor -- ts -->number of rows in tile , ts -->number of columns in tile, ts*ts-->title size (number of elements), N -->number of rows in the whole matrix, N --> number pf columns in the whole matrix
            CHAMELEON_Desc_Create(&CHAM_descC, C, ChamRealDouble, ts, ts, ts * ts, LDC, N, 0,
                    0, LDC, N,1,1);

            CHAMELEON_Desc_Create(&CHAMELEON_descZ, Z, ChamRealDouble, ts, ts,ts*ts, LDZ, NRHS , 0, 0, LDZ, NRHS,1,1);

            CHAMELEON_Desc_Create(&CHAMELEON_descZcpy, Zcpy, ChamRealDouble, ts, ts, ts*ts, LDZ, NRHS, 0, 0,LDZ, NRHS,1,1);



            // Generate Observations Vector (Z) for testing phase
            if(async==0)
            {
                EXAGEOSTAT_ML_GenerateZVec_Tile(obs_dir,obs_timestamp,&r[i*N]);
                START_TIMING(time_opt);
                EXAGEOSTAT_MLE_Tile(1,starting_theta,NULL,NULL);
                //Maximum Likelihood (Optimization Phase using dlib library)
                STOP_TIMING(time_opt);
            }

            else
            {
                EXAGEOSTAT_ML_GenerateZVec_Tile_Async(obs_dir,obs_timestamp,&r[i*N]);
                START_TIMING(time_opt);
                
                nlopt_set_max_objective(opt, EXAGEOSTAT_MLE_Tile_Async, NULL);
                nlopt_optimize(opt, starting_theta, &opt_f);
                
                STOP_TIMING(time_opt);
            }


        }

        printf("Total Number of Iterations=%d\n",iter_count);
        printf("Total Optimization Time= %6.2f\n",time_opt);
        printf("Found Maximum at f(%g, %g, %g) \n", starting_theta[0],
                starting_theta[1],starting_theta[2]);
        int j=0;
        for(j=0;j<3;j++)
            avg_target_theta[j]+=starting_theta[j];

        // for experiments
        total_exec_time=time_opt;
        avg_exec_time_per_iter=avg_exec_time_per_iter/iter_count;
        avg_flops_per_iter=avg_flops_per_iter/iter_count;

        if(zvecs==1)
            write_to_file("results.txt",N,ncores,ts,test,ikernel,based_sys,async,obs_dir,obs_timestamp,total_exec_time,avg_exec_time_per_iter,avg_flops_per_iter);

    }

    printf("\n**************************************************************************************\n");
    printf("Average Maximum using (%d) Z vectors at f(%g, %g, %g) \n", zvecs, avg_target_theta[0]/zvecs, avg_target_theta[1]/zvecs,avg_target_theta[2]/zvecs);
    //Free Memory Allocation
    ML_finalize();
    nlopt_destroy(opt);

    return 0;
}



//Generate Observations Vector (Z) for testing Maximum Likelihood function
int PLASMA_ML_GenerateZVec_Tile (char * obs_dir, int obs_timestamp, int * r) {

    //Initialization
    double flops=0;

    if(strcmp(obs_dir,"")==0)

    {

        //Set C diagonal equals to 1 (simple case)
        PLASMA_dlaset_Tile(PlasmaUpperLower, 0, 1, descC);

        if (verbose == 1)
            fprintf(stderr, "Initializing Co-variance Matrix ...\n");

        //Generate co-variance matrix C
        PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], initial_theta);

        if (verbose == 1)
            fprintf(stderr, "Done ...\n");




        //    int IPV;
        //  PLASMA_dgetrf_Tile(descC,&IPV);
        //	double res=0;
        //      PLASMA_dgecon_Tile(PlasmaOneNorm,descC,PlasmaOneNorm,&res);
        //	printf("res: %f\n",res);
        //	exit(0);



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
        PLASMA_MLE_GenZVec_Tile_Async(descZ, r, sequence, &request[1]);

        PLASMA_dtrmm_Tile(PlasmaLeft,PlasmaLower,PlasmaNoTrans,PlasmaNonUnit,1,descC,descZ);
    }	
    else
    {

        if (verbose == 1)
            fprintf(stderr,
                    "Reading observations from dir ...\n");

        FILE * fp;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        char * pch;
        int str_count=0;
        char* zeroes;
        char integer_string[32];
        int i=0;

        char * obs_dir_name=(char *) malloc(50 * sizeof(char));

        for(i=0;i<103;i++)
        {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i+1));

            if(i<9)
                zeroes="00";
            else if(i>=9&&i<99)
                zeroes="0";
            else
                zeroes="";

            strcat(obs_dir_name,"/ppt.complete.Y");

            strcat(obs_dir_name,zeroes);
            strcat(obs_dir_name,integer_string);



            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count=0;

                pch = strtok (line," ,");
                while (pch != NULL)
                {
                    if(str_count==obs_timestamp)
                        Z[i]=atof(pch);
                    pch = strtok (NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }



    //Copy Z to Zcpy
    PLASMA_dlacpy_Tile(PlasmaUpperLower ,descZ,descZcpy);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");




    return 0;

}



//Generate Observations Vector (Z) for testing Maximum Likelihood function
int PLASMA_ML_GenerateZVec_Tile_Async(char * obs_dir, int obs_timestamp, int * r) {

    //Initialization
    double flops=0;


    if(strcmp(obs_dir,"")==0)

    {

        //		PLASMA_dlaset_Tile_Async(PlasmaUpperLower, 0, 1, descC,sequence, &request[1]);


        if (verbose == 1)
            fprintf(stderr, "Initializing Co-variance Matrix ...\n");


        //Generate co-variance matrix C
        PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], initial_theta);


        //exit(0);
        if (verbose == 1)
            fprintf(stderr, "Done ...\n");

        //PLASMA_Sequence_Wait(sequence);
        //print_matrix("nn",N,N,C,N);
        //exit(0);




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
        PLASMA_MLE_GenZVec_Tile_Async(descZ, r, sequence, &request[1]);

        PLASMA_dtrmm_Tile_Async(PlasmaLeft,PlasmaLower,PlasmaNoTrans,PlasmaNonUnit,1,descC,descZ,sequence, &request[1]);
    }	
    else
    {

        if (verbose == 1)
            fprintf(stderr,
                    "Reading observations from dir ...\n");

        FILE * fp;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        char * pch;
        int str_count=0;
        char* zeroes;
        char integer_string[32];
        int i=0;

        char * obs_dir_name=(char *) malloc(50 * sizeof(char));

        for(i=0;i<103;i++)
        {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i+1));

            if(i<9)
                zeroes="00";
            else if(i>=9&&i<99)
                zeroes="0";
            else
                zeroes="";

            strcat(obs_dir_name,"/ppt.complete.Y");

            strcat(obs_dir_name,zeroes);
            strcat(obs_dir_name,integer_string);



            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count=0;

                pch = strtok (line," ,");
                while (pch != NULL)
                {
                    if(str_count==obs_timestamp)
                        Z[i]=atof(pch);
                    pch = strtok (NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }



    //Copy Z to Zcpy
    PLASMA_dlacpy_Tile_Async(PlasmaUpperLower ,descZ,descZcpy,sequence, &request[1]);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");




    return 0;


}


//Generate Observations Vector (Z) for testing Maximum Likelihood function
int EXAGEOSTAT_ML_GenerateZVec_Tile(char * obs_dir, int obs_timestamp,int * r) {

    //Initialization
    double flops=0;

    if(strcmp(obs_dir,"")==0)
    {

        if (verbose == 1)
            fprintf(stderr, "Initializing Co-variance Matrix ...\n");


        //Generate co-variance matrix C
        EXAGEOSTAT_MLE_GenCovMat_Tile_Async(CHAM_descC, msequence, &mrequest[1], initial_theta);

        if (verbose == 1)
            fprintf(stderr, "Done ...\n");


        if (verbose == 1)
            fprintf(stderr,
                    "Cholesky factorization of Sigma (Generation Phase) ...\n");

        //Cholesky factorization
        int success = EXAGEOSTAT_dpotrf_Tile(ChamLower, CHAM_descC);

        if (success != PLASMA_SUCCESS) {
            printf("Factorization cannot be performed..\n"
                    "The matrix is not positive definite\n\n");
            exit(0);
        }

        flops = flops + FLOPS_DPOTRF(N);
        if (verbose == 1)
            fprintf(stderr, " Done.\n");


        // Generate Z Vector
        EXAGEOSTAT_MLE_GenZVec_Tile_Async(CHAMELEON_descZ, r, msequence, &mrequest[1]);

        EXAGEOSTAT_dtrmm_Tile(ChamLeft,ChamLower,ChamNoTrans,ChamNonUnit,1,CHAM_descC,CHAMELEON_descZ);
    }
    else
    {
        if (verbose == 1)
            fprintf(stderr,
                    "Reading observations from dir ...\n");

        FILE * fp;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        char * pch;
        int str_count=0;
        char* zeroes;
        char integer_string[32];
        int i=0;

        char * obs_dir_name=(char *) malloc(50 * sizeof(char));

        for(i=0;i<103;i++)
        {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i+1));

            if(i<9)
                zeroes="00";
            else if(i>=9&&i<99)
                zeroes="0";
            else
                zeroes="";

            strcat(obs_dir_name,"/ppt.complete.Y");

            strcat(obs_dir_name,zeroes);
            strcat(obs_dir_name,integer_string);



            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count=0;

                pch = strtok (line," ,");
                while (pch != NULL)
                {
                    if(str_count==obs_timestamp)
                        Z[i]=atof(pch);
                    pch = strtok (NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }

    //Copy Z to Zcpy
    EXAGEOSTAT_dlacpy_Tile(ChamUpperLower ,CHAMELEON_descZ,CHAMELEON_descZcpy);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");




    return 0;

}

int EXAGEOSTAT_ML_GenerateZVec_Tile_Async(char * obs_dir, int obs_timestamp, int * r) {

    //Initialization
    double flops=0;

    if(strcmp(obs_dir,"")==0)
    {
        //Set C diagonal equals to 1 (simple case)
        if (verbose == 1)
            fprintf(stderr, "Initializing Co-variance Matrix ...\n");


        //Generate co-variance matrix C
        EXAGEOSTAT_MLE_GenCovMat_Tile_Async(CHAM_descC, msequence, &mrequest[1], initial_theta);


        //exit(0);

        if (verbose == 1)
            fprintf(stderr, "Done ...\n");


        if (verbose == 1)
            fprintf(stderr,
                    "Cholesky factorization of Sigma (Generation Phase) ...\n");

        //Cholesky factorization
        int success = EXAGEOSTAT_dpotrf_Tile_Async(ChamLower, CHAM_descC,msequence,&mrequest[1]);

        if (success != PLASMA_SUCCESS) {
            printf("Factorization cannot be performed..\n"
                    "The matrix is not positive definite\n\n");
            exit(0);
        }

        flops = flops + FLOPS_DPOTRF(N);
        if (verbose == 1)
            fprintf(stderr, " Done.\n");


        // Generate Z Vector
        EXAGEOSTAT_MLE_GenZVec_Tile_Async(CHAMELEON_descZ, r,msequence, &mrequest[1]);

        EXAGEOSTAT_dtrmm_Tile_Async(ChamLeft,ChamLower,ChamNoTrans,ChamNonUnit,1,CHAM_descC,CHAMELEON_descZ,msequence,&mrequest[1]);
    }
    else
    {
        if (verbose == 1)
            fprintf(stderr,
                    "Reading observations from dir ...\n");

        FILE * fp;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        char * pch;
        int str_count=0;
        char* zeroes;
        char integer_string[32];
        int i=0;

        char * obs_dir_name=(char *) malloc(50 * sizeof(char));

        for(i=0;i<103;i++)
        {

            strcpy(obs_dir_name, obs_dir);

            sprintf(integer_string, "%d", (i+1));

            if(i<9)
                zeroes="00";
            else if(i>=9&&i<99)
                zeroes="0";
            else
                zeroes="";

            strcat(obs_dir_name,"/ppt.complete.Y");

            strcat(obs_dir_name,zeroes);
            strcat(obs_dir_name,integer_string);



            fp = fopen(obs_dir_name, "r");
            if (fp == NULL)
                exit(EXIT_FAILURE);

            while ((read = getline(&line, &len, fp)) != -1) {
                str_count=0;

                pch = strtok (line," ,");
                while (pch != NULL)
                {
                    if(str_count==obs_timestamp)
                        Z[i]=atof(pch);
                    pch = strtok (NULL, " ,");
                    str_count++;
                }


            }

        }


        if (line)
            free(line);
    }

    //Copy Z to Zcpy
    EXAGEOSTAT_dlacpy_Tile_Async(ChamUpperLower ,CHAMELEON_descZ,CHAMELEON_descZcpy,msequence,&mrequest[1]);


    if (verbose == 1)
        fprintf(stderr, "Done Z Vector Generation Phase.\n");




    return 0;


}


double PLASMA_MLE_Tile(unsigned n, const double * theta, double * grad, void * my_func_data)
{


    //printf("******************8hi\n");
    //exit(0);
    //Initialization
    double theta_hat,det=1.0,logdet=0.0;
    double matrix_gen_time=0.0;	double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0;
    double flops = 0.0;




    //	PLASMA_dlaset_Tile(PlasmaUpperLower, 0, 1, descC);

    START_TIMING(matrix_gen_time);
    //Generate new co-variance matrix C based on new theta
    PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], theta);
    STOP_TIMING(matrix_gen_time);

    //re-store old Z
    PLASMA_dlacpy_Tile(PlasmaUpperLower ,descZcpy,descZ);


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



    PLASMA_MLE_det_Tile(descC, sequence, &request[1],&det);


    logdet= det==0?0:log(det*det);


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
    PLASMA_dtrsm_Tile	(PlasmaLeft,PlasmaLower,PlasmaNoTrans,PlasmaNonUnit,1,descC,descZ);

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
    printf("------logdet: %2.6f ",logdet);
    printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1], theta[2],	theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    printf(" ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  +logdet_calculate+ time_solve));
    printf("***************************************************\n");


    iter_count++;


    // for experiments
    avg_exec_time_per_iter+=matrix_gen_time+time_facto + logdet_calculate + time_solve;
    avg_flops_per_iter+=flops / 1e9 / (time_facto  + logdet_calculate+time_solve);
    final_theta_hat=theta_hat;


    return theta_hat;
}

double PLASMA_MLE_Tile_Async(unsigned n, const double * theta, double * grad, void * my_func_data) {

    //Initialization
    double theta_hat,det=1.0,logdet=0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0;
    double flops = 0.0;




    //	PLASMA_dlaset_Tile_Async(PlasmaUpperLower, 0, 1, descC,sequence,&request[1]);

    //printf("%f: %f:%f \n",theta[0],theta[1],theta[2]);
    //exit(0);
    START_TIMING(matrix_gen_time);
    //Generate new co-variance matrix C based on new theta
    PLASMA_MLE_GenCovMat_Tile_Async(descC, sequence, &request[1], theta);
    STOP_TIMING(matrix_gen_time);

    //	PLASMA_Sequence_Wait(sequence);

    //printf("%f: %f:%f \n",theta[0],theta[1],theta[2]);
    //print_matrix("test:",N,N,C,N);

    //char  x[12];
    //scanf("%s",x);
    //exit(0);
    //re-store old Z
    PLASMA_dlacpy_Tile_Async(PlasmaUpperLower ,descZcpy,descZ,sequence,&request[1]);



    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = PLASMA_dpotrf_Tile_Async(PlasmaLower, descC,sequence,&request[1]);


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
    PLASMA_dtrsm_Tile_Async	(PlasmaLeft,PlasmaLower,PlasmaNoTrans,PlasmaNonUnit,1,descC,descZ,sequence,&request[1]);

    STOP_TIMING(time_solve);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");

    flops = flops + FLOPS_DTRSM(PlasmaLeft,N, NRHS);

    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);




    PLASMA_Sequence_Wait(sequence);
    PLASMA_MLE_det_Tile(descC, sequence, &request[1],&det);



    logdet= det==0? 0:log(det*det);


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
    printf("------logdet: %2.6f ",logdet);
    printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1], theta[2],	theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    printf(" ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  +logdet_calculate+ time_solve));
    printf("***************************************************\n");


    iter_count++;


    // for experiments
    avg_exec_time_per_iter+=matrix_gen_time+time_facto + logdet_calculate + time_solve;
    avg_flops_per_iter+=flops / 1e9 / (time_facto  + logdet_calculate+time_solve);
    final_theta_hat=theta_hat;



    return theta_hat;
}


double EXAGEOSTAT_MLE_Tile(unsigned n, const double * theta, double * grad, void * my_func_data) {

    //Initialization
    double theta_hat,det=1.0,logdet=0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0;
    double flops = 0.0;
    double test_time;

    START_TIMING(test_time);

    START_TIMING(matrix_gen_time);
    //Generate new co-variance matrix C based on new theta
    EXAGEOSTAT_MLE_GenCovMat_Tile_Async(CHAM_descC, msequence, &mrequest[1], theta);
    STOP_TIMING(matrix_gen_time);


    START_TIMING(test_time);
    //re-store old Z
    EXAGEOSTAT_dlacpy_Tile(ChamUpperLower ,CHAMELEON_descZcpy,CHAMELEON_descZ);
    STOP_TIMING(test_time);

    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = EXAGEOSTAT_dpotrf_Tile(ChamLower, CHAM_descC);

    STOP_TIMING(time_facto);

    if (success != PLASMA_SUCCESS) {
        printf("Factorization cannot be performed..\n"
                "The matrix is not positive definite\n\n");

        exit(0);
    }

    //STOP_TIMING(test_time);


    flops = flops + FLOPS_DPOTRF(N);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");


    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);




    EXAGEOSTAT_MLE_det_Tile(CHAM_descC, msequence, &mrequest[1],&det);


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
    EXAGEOSTAT_dtrsm_Tile(ChamLeft,ChamLower,ChamNoTrans,ChamNonUnit,1,CHAM_descC,CHAMELEON_descZ);

    STOP_TIMING(time_solve);

    flops = flops + FLOPS_DTRSM(ChamLeft,N, NRHS);

    if (verbose == 1)
        fprintf(stderr, "Calculating the log likelihood ...");
    theta_hat = -0.5 * cblas_ddot(N, Z, 1, Z, 1) - 0.5 * logdet
        - (double) (N / 2) * log(2 * PI);

    //STOP_TIMING(test_time);


    if (verbose == 1)
        fprintf(stderr, " Done.\n");

    // Print Iteration Summary
    printf("***************************************************\n");
    printf("------logdet: %2.6f ",logdet);
    printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1], theta[2],	theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    printf(" ---- Test Time: %6.2f\n", test_time);

    printf(" ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop: %6.2f\n", flops / 1e9 );	
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  + time_solve));
    printf("***************************************************\n");


    iter_count++;


    // for experiments
    avg_exec_time_per_iter+=matrix_gen_time+time_facto + logdet_calculate + time_solve;
    avg_flops_per_iter+=flops / 1e9 / (time_facto  + logdet_calculate+time_solve);
    final_theta_hat=theta_hat;

    //to be removed
    //if((flops / 1e9 / (time_facto  + time_solve))<700)
    //{
    //exit(0);
    //}

    return theta_hat;
}

double EXAGEOSTAT_MLE_Tile_Async(unsigned n, const double * theta, double * grad, void * my_func_data) {

    //Initialization
    double theta_hat,det=1.0,logdet=0.0;
    double time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time=0.0;
    double flops = 0.0;

    START_TIMING(matrix_gen_time);
    //Generate new co-variance matrix C based on new theta
    EXAGEOSTAT_MLE_GenCovMat_Tile_Async(CHAM_descC, msequence, &mrequest[1], theta);
    STOP_TIMING(matrix_gen_time);

    //re-store old Z
    EXAGEOSTAT_dlacpy_Tile_Async(ChamUpperLower ,CHAMELEON_descZcpy,CHAMELEON_descZ,msequence,&mrequest[1]);



    //Calculate Cholesky Factorization (C=LL-1)
    if (verbose == 1)
        fprintf(stderr, "Cholesky factorization of Sigma...");

    START_TIMING(time_facto);
    int success = EXAGEOSTAT_dpotrf_Tile_Async(ChamLower, CHAM_descC, msequence,&mrequest[1]);


    STOP_TIMING(time_facto);

    if (success != CHAMELEON_SUCCESS) {
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
    EXAGEOSTAT_dtrsm_Tile_Async	(ChamLeft,ChamLower,ChamNoTrans,ChamNonUnit,1,CHAM_descC,CHAMELEON_descZ,msequence,&mrequest[1]);

    STOP_TIMING(time_solve);
    if (verbose == 1)
        fprintf(stderr, " Done.\n");

    flops = flops + FLOPS_DTRSM(ChamLeft,N, NRHS);

    //Calculate log(|C|) --> log(square(|L|))
    if (verbose == 1)
        fprintf(stderr, "Calculating the log determinant ...");

    START_TIMING(logdet_calculate);




    CHAMELEON_Sequence_Wait(msequence);
    EXAGEOSTAT_MLE_det_Tile(CHAM_descC, msequence, &mrequest[1],&det);



    logdet= det==0? 0:log(det*det);


    STOP_TIMING(logdet_calculate);



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
    printf("------logdet: %2.6f ",logdet);
    printf("------expr1: %2.6f ",(-0.5 * cblas_ddot(N, Z, 1, Z, 1)));
    printf("------expr2: %2.6f ",((double) (N / 2) * log(2 * PI)));
    printf(" ---- Theta1: %2.6f ----  Theta2: %2.6f ---- Theta3: %2.6f ----Theta_hat: %2.6f\n", theta[0], theta[1], theta[2],	theta_hat);
    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    printf(" ---- Total Time: %6.2f\n", matrix_gen_time+time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto  +logdet_calculate+ time_solve));
    printf("***************************************************\n");


    iter_count++;


    // for experiments
    avg_exec_time_per_iter+=matrix_gen_time+time_facto + logdet_calculate + time_solve;
    avg_flops_per_iter+=flops / 1e9 / (time_facto  + logdet_calculate+time_solve);
    final_theta_hat=theta_hat;

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


int GenerateXYLoc(int n, char * locs_file, int based_sys) {


    int i;
    
    if(strcmp(locs_file,"")==0)
    {
        //Uniform random generation of distance matrix (For testing Phase)
        if (verbose == 1)
            fprintf(stderr, "Initializing Locations ...\n");

        for (i = 0; i < n; i++) {
            X[i] = (R - 0.5 + uniform_distribution(MIN_RAND, MAX_RAND))
                * sqrt((double) n);
            Y[i] = (L - 0.5 + uniform_distribution(MIN_RAND, MAX_RAND))
                * sqrt((double) n);
        }
    }
    else
    {

        //exit(0);
        if (verbose == 1)
            fprintf(stderr, "Reading Locations from file ...\n");


        FILE * fp;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        char * pch;
        int str_count=0; 

        fp = fopen(locs_file, "r");
        if (fp == NULL)
            exit(EXIT_FAILURE);

        //avoid header
        read = getline(&line, &len, fp); 
        while ((read = getline(&line, &len, fp)) != -1) {
            str_count=0;

            pch = strtok (line," ,");
            while (pch != NULL)
            {
                if(str_count==1)
                    X[i]=atof(pch);
                else if(str_count==2)
                    Y[i]=atof(pch);
                pch = strtok (NULL, " ,");
                str_count++;
            }

            //	printf("X[i]=: %f, Y[i]=: %f\n", X[i],Y[i]);
            i++;

        }

        fclose(fp);
        if (line)
            free(line);

        printf("\nMatrix Size: %d, \n\n",N);
    }

    return 0;
}

int PLASMA_MLE_GenCovMat_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
        PLASMA_request *request, double *  theta) {

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
            //if( n > m) break;
            tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
            QUARK_CORE_dGenCovMat(quark, &task_flags, tempmm, tempnn, A(m, n),ldam, m * A.mb, n * A.nb, theta);
        }

    }

    return PLASMA_SUCCESS;
}



int PLASMA_MLE_GenZVec_Tile_Async(PLASMA_desc *descA, int * r, PLASMA_sequence *sequence,
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






    //mt is the number of tile rows of the sub-matrix -- nt is the number of tile columns of the sub-matrix -- mb the number of rows in a tile
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam= BLKLDD(A, m);


        QUARK_CORE_dGenZVec(quark, &task_flags, tempmm, 1, A(m,0 ),ldam, m * A.mb, 0, r);
    }

    return PLASMA_SUCCESS;
}


int PLASMA_MLE_det_Tile_Async(PLASMA_desc *descA, PLASMA_sequence *sequence,
        PLASMA_request *request,double * det_acc) {


    //Dynamic scheduler functions
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m, n;
    int ldam;
    int tempnn, tempmm;
    PLASMA_desc A = *descA;

    *det_acc=1;


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
            QUARK_CORE_det(quark, &task_flags, tempmm, tempnn, A(m, n),ldam, m * A.mb, n * A.nb, det_acc);


        }

    }


    QUARK_Barrier(quark);

    return PLASMA_SUCCESS;
}


void QUARK_CORE_dGenCovMat(Quark *quark, Quark_Task_Flags *task_flags, int m,
        int n, double *A, int lda, int m0, int n0, double * theta) {
    /*DAG_CORE_GenCovMat;
      Called by the master thread. Add a new task to the scheduler, providing the data pointers, sizes, and dependency information.
      This function provides the main user interface for the user to write data-dependent algorithms.*/
    //printf("%f, %f, %f, \n", theta[0],theta[1],theta[2]);
    //exit(0);


    QUARK_Insert_Task(quark, CORE_dGenCovMat_quark, task_flags, sizeof(int), &m,
            VALUE, sizeof(int), &n, VALUE, sizeof(double) * lda * n, A,
            OUTPUT | LOCALITY, sizeof(int), &lda, VALUE, sizeof(int), &m0,
            VALUE, sizeof(int), &n0, VALUE, 3*sizeof(double), theta, INPUT, 0);

}

void QUARK_CORE_dGenZVec(Quark *quark, Quark_Task_Flags *task_flags, int m,
        int n, double *A, int lda, int m0, int n0, double * r) {
    /*DAG_CORE_GenCovMat;
      Called by the master thread. Add a new task to the scheduler, providing the data pointers, sizes, and dependency information.
      This function provides the main user interface for the user to write data-dependent algorithms.*/
    QUARK_Insert_Task(quark, CORE_dGenZVec_quark, task_flags, sizeof(int), &m,
            VALUE, sizeof(int), &n, VALUE, sizeof(double) * lda * n, A,
            OUTPUT | LOCALITY, sizeof(int), &lda, VALUE, sizeof(int), &m0,
            VALUE, sizeof(int), &n0, VALUE,sizeof(double)*lda*n,r,INPUT, 0);


}



void QUARK_CORE_det(Quark *quark, Quark_Task_Flags *task_flags, int m,
        int n, double *A, int lda, int m0, int n0, double * determinant ) {
    /*DAG_CORE_GenCovMat;
      Called by the master thread. Add a new task to the scheduler, providing the data pointers, sizes, and dependency information.
      This function provides the main user interface for the user to write data-dependent algorithms.*/
    QUARK_Insert_Task(quark, CORE_det_quark, task_flags, sizeof(int), &m,
            VALUE, sizeof(int), &n, VALUE, sizeof(double) * lda * n, A,
            INPUT, sizeof(int), &lda, VALUE, sizeof(int), &m0,
            VALUE, sizeof(int), &n0, VALUE,  sizeof(double), determinant, INOUT, 0);

}

void CORE_dGenCovMat_quark(Quark *quark) {

    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double * theta;


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
    double * r;

    quark_unpack_args_7(quark, m, n, A, lda, m0, n0,r);
    veccov_comp_Tile(A, m, n, m0, n0,r);

}


void  CORE_det_quark(Quark *quark) {

    int m;
    int n;
    double *A;
    int lda;
    int m0;
    int n0;
    double *determinant;

    quark_unpack_args_7(quark, m, n, A, lda, m0, n0, determinant);
    det_comp_Tile(A, m, n, m0, n0,determinant);
}

void matcov_comp_Tile(double * A, int m, int n, int m0, int n0, double * theta) {



    int i, j;
    int x = m0;
    int y = n0;

    double dist=0.0;
    double expr=0.0;
    double con=0.0;



    con=pow(2,(theta[2]-1)) * tgamma(theta[2]);

    con=1.0/con;
    con=theta[0]*con;

    for (i = 0; i < m; i++) {
        y = n0;

        for (j = 0; j < n; j++) {

            expr=CalculateDistance(X[x], Y[x], X[y], Y[y])/theta[1];
            expr= expr==0? 1e-10:expr;
            A[i + j * m]=con*pow(expr,theta[2])*gsl_sf_bessel_Knu(theta[2],expr); // Matern Function

            //if(r1!=r2)
            //printf("%f - %f - %f - %f - %f - %f - %f - %f - %f\n",A[i + j * m],X[x],Y[x],X[y],Y[y],theta[1],con,expr,gsl_sf_bessel_Knu(theta[2],expr));

            y++;
        }
        x++;
    }

    /*

       for (i = 0; i < m; i++) {
       y = n0;
       for (j = 0; j < n; j++) {
       A[i + j * m] =
       A[i + j * m] == 1 ?
1 :
exp(-CalculateDistance(X[x], Y[x], X[y], Y[y]) / theta[1]);

y++;
}
x++;
}
*/
}





void det_comp_Tile(double * A, int m, int n, int m0, int n0,double * res) {

    //*res=1;
    //printf("yes\n");
    int i, j;
    int x = m0;
    int y = n0;

    for (i = 0; i < m; i++) {
        y = n0;
        for (j = 0; j < n; j++) {
            if(x==y)
            {
                *res*=A[i + j * m];
                if(*res==0)
                    return;
                //printf("A[i+j*m: %f \n",A[i+j*m]);
            }

            y++;
        }
        x++;
    }

    //printf("*res= %f\n", *res);
}




void veccov_comp_Tile(double * Z, int m, int n, int m0, int n0, double *r) {

    int i=0;

    for(i=0;i<m;i++)
        Z[i]=r[m0+i];
}

void print_matrix(char* desc, int m, int n, double* a, int lda) {
    int i, j;
    printf("\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf(" %6.4f", a[i + lda * j]);
        printf("\n");
    }
}

//*******************************************************************************(1)
void CORE_dGenCovMat_starpu(void *buffers[],void *cl_arg){
    int m,n,lda,m0,n0;
    double * theta;
    double *A;

    starpu_codelet_unpack_args(cl_arg, &m,&n,&A,&lda,&m0,&n0,&theta);
    matcov_comp_Tile(A, m, n, m0, n0, theta);
}

static struct starpu_codelet cl_dGenCovMat =
{
    .where = STARPU_CPU,
    .cpu_funcs = {CORE_dGenCovMat_starpu},
    .nbuffers = 0
};


int EXAGEOSTAT_MLE_GenCovMat_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence, RUNTIME_request_t  *request, double * theta) {

    CHAM_context_t *morse;
    RUNTIME_option_t options;
    morse = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m, n, m0, n0;
    int ldam;
    int tempmm,tempnn;
    CHAM_desc_t A = *descA;

    struct starpu_codelet *cl=&cl_dGenCovMat;




    for(m=0; m<A.mt;m++)
    {
        tempmm=m==A.mt -1? A.m- m* A.mb: A.mb;
        ldam= mBLKLDD(descA, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt -1? A.n - n * A.nb: A.nb;

            double * data=(double*)morse_getaddr_ccrb(descA,m,n);

            m0= m * A.mb;
            n0= n * A.nb;
            //      printf("m0=%d, ",m0);
            //	printf("n0=%d, \n",n0);

            starpu_insert_task(cl,
                    STARPU_VALUE, &tempmm,  sizeof(int),
                    STARPU_VALUE, &tempnn, sizeof(int),
                    STARPU_VALUE, &data, sizeof(double *),
                    STARPU_VALUE, &ldam,   sizeof(int),
                    STARPU_VALUE, &m0,   sizeof(int),
                    STARPU_VALUE, &n0,   sizeof(int),
                    STARPU_VALUE, &theta,   sizeof(double *),
                    0);
        }
    }

    RUNTIME_options_finalize(&options, morse);
    CHAMELEON_TASK_dataflush_all();
    return CHAMELEON_SUCCESS;
}

//****************************************************************************(2)
void CORE_GenZVec_starpu(void *buffers[],void *cl_arg){
    int m;
    int n;
    double * A;
    int lda;
    int m0;
    int n0;
    double * r;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &A, &lda, & m0, &n0, &r );
    veccov_comp_Tile(A, m, n, m0, n0, r);
}

static struct starpu_codelet cl_GenZVec =
{
    .where = STARPU_CPU,
    .cpu_funcs = {CORE_GenZVec_starpu},
    .nbuffers = 0
};


int EXAGEOSTAT_MLE_GenZVec_Tile_Async(CHAM_desc_t *descA, int * r,RUNTIME_sequence_t *sequence, RUNTIME_request_t  *request) {

    CHAM_context_t *morse;
    RUNTIME_option_t options;
    morse = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m,m0,n0;
    int ldam;
    int tempmm, tempnn;
    CHAM_desc_t A = *descA;



    struct starpu_codelet *cl=&cl_GenZVec;


    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m - m * A.mb : A.mb;
        ldam = mBLKLDD(descA, m);

        tempnn=1;

        double *data=(double*)morse_getaddr_ccrb(descA,m,0);

        m0= m * A.mb;
        n0= 0;

        starpu_insert_task(cl,
                STARPU_VALUE, &tempmm,  sizeof(int),
                STARPU_VALUE, &tempnn, sizeof(int),
                STARPU_VALUE, &data,   sizeof(double*),
                STARPU_VALUE, &ldam,   sizeof(int),
                STARPU_VALUE, &m0,   sizeof(int),
                STARPU_VALUE, &n0,   sizeof(int),
                STARPU_VALUE, &r,   sizeof(double),
                0);

    }

    RUNTIME_options_finalize(&options, morse);
    CHAMELEON_TASK_dataflush_all();
    return CHAMELEON_SUCCESS;
}


//*****************************************************************************(3)
void CORE_det_starpu(void *buffers[],void *cl_arg){
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


int EXAGEOSTAT_MLE_det_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence, RUNTIME_request_t  *request, double * det_acc) {

    CHAM_context_t *morse;
    RUNTIME_option_t options;
    morse = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, morse, sequence, request);

    int m,n,m0,n0;
    int ldam;
    int tempmm, tempnn;
    CHAM_desc_t A = *descA;
    *det_acc=1.0;
    //double det_val=*det_acc;
    struct starpu_codelet *cl=&cl_det;

    for(m=0; m<A.mt;m++)
    {
        tempmm=m==A.mt -1? A.m- m* A.mb: A.mb;
        ldam= mBLKLDD(descA, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt -1? A.n - n * A.nb: A.nb;

            double * data=(double*)morse_getaddr_ccrb(descA,m,n);

            m0= m * A.mb;
            n0= n * A.nb;

            starpu_insert_task(cl,
                    STARPU_VALUE, &tempmm,  sizeof(int),
                    STARPU_VALUE, &tempnn, sizeof(int),
                    STARPU_VALUE, &data,   sizeof(double*),
                    STARPU_VALUE, &ldam,   sizeof(int),
                    STARPU_VALUE, &m0,   sizeof(int),
                    STARPU_VALUE, &n0,   sizeof(int),
                    STARPU_VALUE, &det_acc,   sizeof(double*),
                    0);
            //printf("det_val: %f\n",det_val);
        }
    }
    RUNTIME_barrier(morse);
    RUNTIME_options_finalize(&options, morse);
    CHAMELEON_TASK_dataflush_all();
    return CHAMELEON_SUCCESS;
}

int countlines(char *filename)
{
    // count the number of lines in the file called filename
    FILE *fp = fopen(filename,"r");
    int ch=0;
    int lines=0;

    if (fp == NULL)
    {
        printf("cannot open locations file\n");
        return 0;
    }


    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
            lines++;
    }


    fclose(fp);

    //Excluding header line
    return (lines-1);
}

static double bessy0( double x )
    /*------------------------------------------------------------*/
    /* PURPOSE: Evaluate Bessel function of second kind and order */
    /*          0 at input x.                                     */
    /*------------------------------------------------------------*/
{
    double z;
    double xx,y,ans,ans1,ans2;

    if (x < 8.0) {
        y=x*x;
        ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
                    +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
        ans2=40076544269.0+y*(745249964.8+y*(7189466.438
                    +y*(47447.26470+y*(226.1030244+y*1.0))));
        ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
    } else {
        z=8.0/x;
        y=z*z;
        xx=x-0.785398164;
        ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
                    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
        ans2 = -0.1562499995e-1+y*(0.1430488765e-3
                +y*(-0.6911147651e-5+y*(0.7621095161e-6
                        +y*(-0.934945152e-7))));
        ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
    }
    return ans;
}



static double bessy1( double x )
    /*------------------------------------------------------------*/
    /* PURPOSE: Evaluate Bessel function of second kind and order */
    /*          1 at input x.                                     */
    /*------------------------------------------------------------*/
{
    double z;
    double xx,y,ans,ans1,ans2;

    if (x < 8.0) {
        y=x*x;
        ans1=x*(-0.4900604943e13+y*(0.1275274390e13
                    +y*(-0.5153438139e11+y*(0.7349264551e9
                            +y*(-0.4237922726e7+y*0.8511937935e4)))));
        ans2=0.2499580570e14+y*(0.4244419664e12
                +y*(0.3733650367e10+y*(0.2245904002e8
                        +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
        ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
    } else {
        z=8.0/x;
        y=z*z;
        xx=x-2.356194491;
        ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
                    +y*(0.2457520174e-5+y*(-0.240337019e-6))));
        ans2=0.04687499995+y*(-0.2002690873e-3
                +y*(0.8449199096e-5+y*(-0.88228987e-6
                        +y*0.105787412e-6)));
        ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
    }
    return ans;
}



/*

   Return the Bessel function of second kind and
   of integer order, for input value x.
   n        Integer order of Bessel function.
   x        Double at which the function is evaluated.

*/


double bessy( int n, double x )
    /*------------------------------------------------------------*/
    /* PURPOSE: Evaluate Bessel function of second kind and order */
    /*          n for input x. (n >= 0)                           */
    /* Note that for x == 0 the functions bessy and bessk are not */
    /* defined and a blank is returned.                           */
    /*------------------------------------------------------------*/
{
    int j;
    double by,bym,byp,tox;


    if (n < 0 || x == 0.0)
    {
        double   dblank;
        // setdblank_c( &dblank );
        return( dblank );
    }
    if (n == 0)
        return( bessy0(x) );
    if (n == 1)
        return( bessy1(x) );

    tox=2.0/x;
    by=bessy1(x);
    bym=bessy0(x);
    for (j=1;j<n;j++) {
        byp=j*tox*by-bym;
        bym=by;
        by=byp;
    }
    return by;
}

static double bessj0( double x )
    /*------------------------------------------------------------*/
    /* PURPOSE: Evaluate Bessel function of first kind and order  */
    /*          0 at input x                                      */
    /*------------------------------------------------------------*/
{
    double ax,z;
    double xx,y,ans,ans1,ans2;

    if ((ax=fabs(x)) < 8.0) {
        y=x*x;
        ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
                    +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
        ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
                    +y*(59272.64853+y*(267.8532712+y*1.0))));
        ans=ans1/ans2;
    } else {
        z=8.0/ax;
        y=z*z;
        xx=ax-0.785398164;
        ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
                    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
        ans2 = -0.1562499995e-1+y*(0.1430488765e-3
                +y*(-0.6911147651e-5+y*(0.7621095161e-6
                        -y*0.934935152e-7)));
        ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    }
    return ans;
}



static double bessj1( double x )
    /*------------------------------------------------------------*/
    /* PURPOSE: Evaluate Bessel function of first kind and order  */
    /*          1 at input x                                      */
    /*------------------------------------------------------------*/
{
    double ax,z;
    double xx,y,ans,ans1,ans2;

    if ((ax=fabs(x)) < 8.0) {
        y=x*x;
        ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
                        +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
        ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
                    +y*(99447.43394+y*(376.9991397+y*1.0))));
        ans=ans1/ans2;
    } else {
        z=8.0/ax;
        y=z*z;
        xx=ax-2.356194491;
        ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
                    +y*(0.2457520174e-5+y*(-0.240337019e-6))));
        ans2=0.04687499995+y*(-0.2002690873e-3
                +y*(0.8449199096e-5+y*(-0.88228987e-6
                        +y*0.105787412e-6)));
        ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
        if (x < 0.0) ans = -ans;
    }
    return ans;
}

/*int PLASMA_MLE_det_Tile(PLASMA_desc *descA, PLASMA_sequence *sequence,
  PLASMA_request *request,double * det_acc) {


  int m,m0,n0,n;
  int ldam;
  int tempnn, tempmm;
  PLASMA_desc A = *descA;

 *det_acc=1;



//mt is the number of tile rows of the sub-matrix -- nt is the number of tile columns of the sub-matrix -- mb the number of rows in a tile
for (m = 0; m < A.mt; m++) {
tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
ldam = BLKLDD(A, m);

for (n = 0; n < A.nt; n++) {
//generate the Lower and diagonal tiles if symmetric
// if(part == 1 && n > m) break;
tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;

m0=m*A.mb;
n0=n*A.nb;

det_comp_Tile(A(m,n), tempmm, tempnn, m0, n0,det_acc);



}

}


PLASMA_Sequence_Wait();


return PLASMA_SUCCESS;
}
*/





int PLASMA_MLE_det_Tile(PLASMA_desc *descA, PLASMA_sequence *sequence,
        PLASMA_request *request,double * det_acc) {


    //Dynamic scheduler functions
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m, n, m0, n0;
    int ldam;
    int tempnn, tempmm;
    PLASMA_desc A = *descA;

    *det_acc=1;


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

            m0=m*A.mb;
            n0=n*A.nb;

            det_comp_Tile(A(m,n), tempmm, tempnn, m0, n0,det_acc);


        }

    }


    QUARK_Barrier(quark);

    return PLASMA_SUCCESS;
}


int EXAGEOSTAT_MLE_det_Tile(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence, RUNTIME_request_t  *request, double * det_acc) {

    CHAM_context_t *morse;
    RUNTIME_option_t options;
    morse = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, morse, sequence, request);

    int m,n,m0,n0;
    int ldam;
    int tempmm, tempnn;
    CHAM_desc_t A = *descA;
    *det_acc=1.0;
    //double det_val=*det_acc;
    struct starpu_codelet *cl=&cl_det;

    for(m=0; m<A.mt;m++)
    {
        tempmm=m==A.mt -1? A.m- m* A.mb: A.mb;
        ldam= mBLKLDD(descA, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt -1? A.n - n * A.nb: A.nb;

            double * data=(double*)morse_getaddr_ccrb(descA,m,n);

            m0= m * A.mb;
            n0= n * A.nb;

            det_comp_Tile(data, tempmm, tempnn, m0, n0,det_acc);
            //printf("det_val: %f\n",det_val);
        }
    }
    RUNTIME_barrier(morse);

    RUNTIME_options_finalize(&options, morse);
    CHAMELEON_TASK_dataflush_all();
    return CHAMELEON_SUCCESS;
}



void write_to_file(char * path, int matrix_size,int ncores,int tile_size, int test, char * ikernel, int based_system, int async, char *obs_dir,int obs_timestamp,double total_exec_time,double avg_exec_time_per_iter, double avg_flops_per_iter)
{
    FILE *pFile;

    pFile=fopen(path,"a");

    if(pFile==NULL) {


        //	    fopen(path, "w+");
        //          fprintf(pFile, "ts-");
        //	    fprintf(pFile, "kernel-");
        //	    fprintf(pFile, "BASED_SYS-");
        //          fprintf(pFile, "SYNC-");

        //	    if(test==1)
        //      	fprintf(pFile, "N\t\t");
        //	    else
        //	      {
        //      	fprintf(pFile, "obs_dir-");
        //	        fprintf(pFile, "time_stamp \t\t");
        //	      }
        //	    fprintf(pFile, "ncores\t");
        //	    fprintf(pFile, "ET\t");
        //	    fprintf(pFile, "ET_Iter\t");
        //	    fprintf(pFile, "FLOPS_Iter\n");

        printf("Cannot access the results path\n");
        exit(0);
    }

    fprintf(pFile, "%d\t", N);
    fprintf(pFile, "%d\t", ncores);
    fprintf(pFile, "%f\t", total_exec_time);
    fprintf(pFile, "%f\t", avg_exec_time_per_iter);
    fprintf(pFile, "%f\t\t", avg_flops_per_iter);



    fprintf(pFile, "%d-", tile_size);
    fprintf(pFile, "%s-", ikernel);

    if(based_system==0)
        fprintf(pFile, "PLASMA-");
    else
        fprintf(pFile, "CHAMELEON-");


    if(async=0)
        fprintf(pFile, "SYNC-");
    else
        fprintf(pFile, "ASYNC-");



    if(test==1)
        fprintf(pFile, "%d-", matrix_size);
    else
    {
        fprintf(pFile, "%s-", obs_dir);
        fprintf(pFile, "%d-", obs_timestamp);
    }

    fprintf(pFile, "%f \n", final_theta_hat);



    fclose(pFile);

}


void theta_parser(double * theta_vec,char * kern)
{

    int i=0;

    if(!strcmp(kern,""))
    {

        for(i=0;i<3;i++)
            theta_vec[i]=-1;

    }


    char * token;
    while( (token = strsep(&kern,":")) != NULL )
    {
        if (strcmp(token,"?"))

            theta_vec[i]=strtod(token,NULL);

        else
            theta_vec[i]=-1;
        i++;
    }


}
