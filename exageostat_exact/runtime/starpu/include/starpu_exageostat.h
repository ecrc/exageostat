/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file starpu_exageostat.h
 *
 * StarPU codelets functions header file.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-19
 *
 **/
#ifndef _EXAGEOSTATCODELETS_H_
#define _EXAGEOSTATCODELETS_H_
#include "../../../core/include/exageostatcore.h"
#if defined(EXAGEOSTAT_USE_CUDA)
#include "../../../../cudacore/include/exageostatcudacore.h"
#endif
#include "morse.h"
#include "../../../../include/morse_starpu.h"
#include "../../../../include/context.h"

//#include "../../../../include/coreblas_ds.h"

//#include "chameleon_starpu.h"
#if defined(CHAMELEON_USE_MPI)
#undef STARPU_REDUX
#define CUBLAS_STREAM_PARAM cublasHandle_t handle
#define starpu_insert_task starpu_mpi_insert_task

#define starpu_mpi_codelet(_codelet_) MPI_COMM_WORLD, _codelet_
#else

#define starpu_mpi_codelet(_codelet_) _codelet_
#endif

#define EXAGEOSTAT_RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)EXAGEOSTAT_data_getaddr( desc, type, m, n ) )

#define mBLKLDD(A, k) A->get_blkldd( A,k )

void *EXAGEOSTAT_data_getaddr( const MORSE_desc_t *A, int type, int m, int n );

int MORSE_MLE_dcmg_Tile_Async(MORSE_enum uplo, MORSE_desc_t *descA, location *l1,
        location *l2, location *lm,  double *theta,
        char *dm, char *kernel_fun,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

int MORSE_MLE_scmg_Tile_Async(MORSE_enum uplo, MORSE_desc_t *descA, location *l1,
        location *l2, location *lm,  double *theta,
        char *dm, char *kernel_fun,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

int MORSE_MLE_sdcmg_Tile_Async(MORSE_enum uplo, MORSE_desc_t *descA, location *l1,
        location *l2, location *lm,  double *theta,
        char *dm, char *kernel_fun,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);


int MORSE_MLE_dzcpy_Tile_Async(MORSE_desc_t *descA, double *r,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

int MORSE_MLE_szcpy_Tile_Async(MORSE_desc_t *descA, float  *r,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

int MORSE_MLE_dmdet_Tile_Async(MORSE_desc_t *descA, 
        MORSE_sequence_t *sequence, MORSE_request_t  *request,
        MORSE_desc_t *descdet);

int MORSE_MLE_smdet_Tile_Async(MORSE_desc_t *descA,
        MORSE_sequence_t *sequence, MORSE_request_t  *request,
        MORSE_desc_t *descdet);

int MORSE_MLE_ddotp_Async(MORSE_enum transA, MORSE_enum transB,
        MORSE_desc_t *descA, MORSE_desc_t *descB,
        MORSE_desc_t *descC, MORSE_enum indexC,
        MORSE_sequence_t *sequence, MORSE_request_t *request);

int MORSE_MLE_dmse_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss,
        MORSE_desc_t *descserror,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

int MORSE_MLE_smse_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss,
        MORSE_desc_t *descserror, 
        MORSE_sequence_t *sequence, MORSE_request_t *request);

int MORSE_MLE_dgemv_Tile_Async(MORSE_desc_t *descA, MORSE_desc_t *descZ,
        MORSE_desc_t *descZout,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

int MORSE_MLE_sgemv_Tile_Async(MORSE_desc_t *descA, MORSE_desc_t *descZ,
        MORSE_desc_t *descZout,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

void MORSE_TASK_sdconv(const MORSE_option_t *options,
        int m, int n, int nb,
        const MORSE_desc_t *A, int Am, int An, int lda,
        const MORSE_desc_t *B, int Bm, int Bn, int ldb);

void MORSE_TASK_dsconv(const MORSE_option_t *options,
        int m, int n, int nb,
        const MORSE_desc_t *A, int Am, int An, int lda,
        const MORSE_desc_t *B, int Bm, int Bn, int ldb);



/*void MORSE_TASK_slag2d(MORSE_option_t *options,
  int m, int n, int nb,
  MORSE_desc_t *A, int Am, int An, int lda,
  MORSE_desc_t *B, int Bm, int Bn, int ldb);
  void MORSE_TASK_dlag2s(MORSE_option_t *options,
  int m, int n, int nb,
  MORSE_desc_t *A, int Am, int An, int lda,
  MORSE_desc_t *B, int Bm, int Bn, int ldb);
  */

/*int MORSE_MLE_dprint_Tile_Async(MORSE_desc_t *descA,
  MORSE_sequence_t *sequence, MORSE_request_t  *request);

  int MORSE_MLE_sprint_Tile_Async(MORSE_desc_t *descA,
  MORSE_sequence_t *sequence, MORSE_request_t  *request);

  int MORSE_MLE_sdprint_Tile_Async(MORSE_desc_t *descA,
  MORSE_sequence_t *sequence, MORSE_request_t *request,
  int diag_thick);
  */

int MORSE_MLE_dmse_bivariate_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss,
        MORSE_desc_t *descserror1,  MORSE_desc_t *descserror2,  MORSE_desc_t *descserror,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

void MORSE_TASK_sexageostat_gemm(const MORSE_option_t *options,
        MORSE_enum transA, int transB,
        int m, int n, int k, int nb,
        float alpha, const MORSE_desc_t *A, int Am, int An, int lda,
        const MORSE_desc_t *B, int Bm, int Bn, int ldb,
        float beta,  const MORSE_desc_t *C, int Cm, int Cn, int ldc);

void MORSE_TASK_sexageostat_trsm(const MORSE_option_t *options,
        MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
        int m, int n, int nb,
        float alpha, const MORSE_desc_t *A, int Am, int An, int lda,
        const MORSE_desc_t *B, int Bm, int Bn, int ldb);


void MORSE_MLE_sdmat_register_Tile_Async(MORSE_enum uplo, MORSE_desc_t *descA,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);

#endif /* _EXAGEOSTATCODELETS_H_ */
