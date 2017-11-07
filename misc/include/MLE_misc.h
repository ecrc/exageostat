/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_misc.h
 *
 * Header file of  auxiliary functions that are needed by ExaGeoStat.
 *
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#ifndef _MLE_MISC_H_
#define _MLE_MISC_H_
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <nlopt.h>
#include <math.h>
#include <morse.h>
#include <starpu.h>
#include <cblas.h>
#include "../../include/flops.h"
#include <starpu_profiling.h>
#if defined(CHAMELEON_USE_MPI)
#include <starpu_mpi.h>
#else
#endif
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>
#endif
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include "../../include/morse_starpu.h"
/** ****************************************************************************
* PI value
**/
#define PI (3.141592653589793)
/** ****************************************************************************
* The radius of the  eartch (used by Great Circle Distance (GCD)
**/
#define earthRadiusKm 6371.0
/** ****************************************************************************
* Start timing macro
**/
#define START_TIMING(_t)     _t =- cWtime();
/** ****************************************************************************
* Stop timing macro
**/
#define STOP_TIMING(_t)      _t += cWtime();


/** ****************************************************************************
 *  EXAGEOSTAT location uniquely identifies two dimensions struct with two vectors (x and y).
 **/
typedef struct {
double *x;				///< x dimension vector.
double *y;				///< y dimension vector.
} location;


/** ****************************************************************************
 *  EXAGEOSTAT accuracy struct  uniquely identifies three different parameters for accuracy
 *  measurments accuracyDenseAppDiff, normDenseAppDiff, and normA.
 **/
typedef struct{
double accuracyDenseAppDiff; 		///< Accuracy different between dense format and approx format.
double normDenseAppDiff;		///< Norm difference between dense format and approx format.
double normA;   			///< normA.
} acc_struct;


/** ****************************************************************************
 *  EXAGEOSTAT MLE  struct  uniquely identifies different variables that are needed
 *  by EXAGEOSTAT library to ease arguments pass function
 **/
typedef struct
{
        char *computation;		///< exact or approx computation.
        int async;			///< running mode: synchronous or asynchronous.
        int iter_count;			///< number of iterations to converge.
        location  l1;			///< 2D locations for the first dataset.
        location  l2;			///< 2D locations for the second dataset.
        void *descC;			///< Covariance matrix C descriptor.
        void *descZ;			///< Measurements Z descriptor.
	double *Adense;			///< Dense matrix descriptor in the case of approximation mode - accuracy check.	
	double *Adense2;		///< Dense matrix descriptor2 in the case of approximation mode - accuracy check.
        void *descZcpy;			///< A copy of Measurements Z descriptor.
        void *descdet;			///< Determinant descriptor.
        void *descproduct;		///< Dot product descriptor.
        void *descZmiss;		///< Missing measurements descriptor.
        void *descC12;			///< Covariance Matrix C12 descriptor.
        void *descC22;			///< Covariance Matrix C22 descriptor.
        void *descZactual;		///< Actual Measurements Z descriptor.
        void *descZobs;			///< observed Measurements Z descriptor.
        void *descmse;         		///< Mean Square Error (MSE) descriptor.
        void *sequence;			///< MORSE sequence.
        void *request;			///< MORSE request
        int verbose;			///< Verbose indicator.
	int check;			///< Check indicator -- approximation mode.
        int log;			///< Log files generation indicator, 0-->no, 1-->yes
        double avg_exec_time_per_iter;	///< Avergae execution time per iteration (only used in verbose mode).
        double total_exec_time;		///< Total execution time (only used in verbose mode).
        double avg_flops_per_iter;	///< Avergae flops per iteration (only used in verbose mode).
        double final_loglik;		///< Final log likelihood value.
        char *locsFPath;		///< Locations file path -- in the case of real dataset (real mode). 
        char *obsFPath;			///< Observations file path --  in the case of real dataset (real mode).
        char *actualZFPath;		///< Actual observations file path -- in the case of prediction.
        char *actualZLocFPath;		///< Actial locations file path -- in the case of prediction.
        double det;			///< determinant value.
        double dotp;			///< dot product value.
        double mserror;			///< Mean Square Error (MSE) value.
        char *dm;			///< Distance metric to be used ed->Euclidian Distance -- gcd->Greate Circle Distance.
        char *nFileLog;			///< log file name (only used if log -->1).
        FILE *pFileLog;			///< log file path (only used if log -->1).
        int hicma_maxrank;		///< Max Rank in the case of LR-HiCMA approx
        //hicma desc_c			 
        void *hicma_descCD;		///< HiCMA descCD descriptor.
        void *hicma_descCUV;		///< HiCMA descCUV descriptor.
        void *hicma_descCrk;		///< HiCMA descCrk descriptor.
        //hicma descC12			///<
        void *hicma_descC12D;		///< HiCMA descCD descriptor.	
        void *hicma_descC12UV;		///< HiCMA descCUV descriptor.
        void *hicma_descC12rk;		///< HiCMA descCrk descriptor.
        //hicma descC12			///<
        void *hicma_descC22D;		///< HiCMA descCD descriptor.	
        void *hicma_descC22UV;		///< HiCMA descCUV descriptor.
	void *hicma_descC22rk;		///< HiCMA descCrk descriptor.
        double hicma_acc;		///< Accuracy in the case of LR-HiCMA approx

} MLE_data;

/** ****************************************************************************
 * Verbose Macro.
 **/
#define VERBOSE(str)	\
        if (data->verbose == 1 && MORSE_My_Mpi_Rank() == 0){	\
        fprintf(stdout, "%s", str);	\
        fflush(stdout);\
    }

/** ****************************************************************************
 *  Success Macro.
 **/
#define SUCCESS(success, str) \
        if (success != MORSE_SUCCESS){ \
        fprintf(stdout, "%s", str);\
        fflush(stdout);\
        exit(EXIT_FAILURE);\
    }




double uniform_distribution(double rangeLow, double rangeHigh);
int GenerateXYLoc(int n, char * locs_file, location * locations);
void print_matrix(char* desc, int m, int n, double* a, int lda);
int countlines(char *filename);
void write_to_file(char * path, int matrix_size,int ncores,int tile_size, int test, char * ikernel, char *computation, int async, char *obsFPath,double total_exec_time,double avg_exec_time_per_iter, double avg_flops_per_iter , int p_grid, int q_grid, double final_loglik, int n);
void theta_parser2(double * theta_vec,char * kern);
void write_vectors(double * zvec, MLE_data * data, int n);
void write_to_thetafile(char * path, double theta0,double theta1, double theta2, double loglik, int n);
void readObsFile(char *obsfile, int n, double * streamdata);
void shuffle(double *array, location * locations, size_t n);
void theta_parser(double *initial_theta, double *target_theta, double *starting_theta, char *ikernel, char *kernel, double *lb, double *up, int test);
void init_optimizer( nlopt_opt * opt, double *lb, double *up, double tol);
void print_summary(int test, int N, int ncores, int gpus, int ts, char *computation, int zvecs, int p_grid, int q_grid);
void print_result(MLE_data *data, double *starting_theta, int N, int zvecs, int ncores, int ts, int test, char * ikernel, char *computation, int p_grid, int q_grid, double final_loglik);
double cWtime(void);
void readlocfile(char* loc_file, int n,  location* l1);
void write_prediction_result(char * path, int matrix_size,int no_missing, double MSE, double solve_time, double flops);
int doesFileExist(const char *filename);
void init_log (MLE_data * data);
void finalize_log (MLE_data * data);




#endif
