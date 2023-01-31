/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
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
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2020-03-04
 *
 **/
#ifndef _EXAGEOSTATCODELETS_H_
#define _EXAGEOSTATCODELETS_H_

#include "../../../cpu_core/include/exageostatcore.h"

#if defined(EXAGEOSTAT_USE_CUDA)
#endif

#include "chameleon.h"
#include "../../../../include/chameleon_starpu.h"
#include "../../../../include/context.h"

#if defined(CHAMELEON_USE_MPI)
#else

#define starpu_mpi_codelet(_codelet_) _codelet_
#endif

#define EXAGEOSTAT_RTBLKADDR(desc, type, m, n) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )
#ifdef EXAGEOSTAT_USE_HICMA

#include <hicma_runtime.h>

#define EXAGEOSTAT_HICMA_RTBLKADDR(desc, type, m, n) ( (starpu_data_handle_t)HICMA_RUNTIME_data_getaddr( desc, m, n ) )

#endif

#define BLKLDD(A, k) A->get_blkldd( A,k )
#define mBLKLDD(A, k) A->get_blkldd( A,k )

void *EXAGEOSTAT_data_getaddr(const CHAM_desc_t *A, int type, int m, int n);

int EXAGEOSTAT_MLE_dcmg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                  location *l2, location *lm, double* theta,
                                  char* dm, char* kernel_fun,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_scmg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                  location *l2, location *lm, double* theta,
                                  char* dm, char* kernel_fun,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_sdcmg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                   location *l2, location *lm, double* theta,
                                   char* dm, char* kernel_fun,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);


int EXAGEOSTAT_MLE_dzcpy_Tile_Async(CHAM_desc_t *descA, double* r,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_szcpy_Tile_Async(CHAM_desc_t *descA, double* r,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dmdet_Tile_Async(CHAM_desc_t *descA,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request,
                                   CHAM_desc_t *descdet);

int EXAGEOSTAT_MLE_smdet_Tile_Async(CHAM_desc_t *descA,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request,
                                   CHAM_desc_t *descdet);

int EXAGEOSTAT_MLE_ddotp_Async(CHAM_enum transA, CHAM_enum transB,
                              CHAM_desc_t *descA, CHAM_desc_t *descB,
                              CHAM_desc_t *descC, CHAM_enum index, CHAM_enum increment,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dmse_Tile_Async(CHAM_desc_t *descZpre, CHAM_desc_t *descZmiss,
                                  CHAM_desc_t *descserror,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_smse_Tile_Async(CHAM_desc_t *descZpre, CHAM_desc_t *descZmiss,
                                  CHAM_desc_t *descserror,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dgemv_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descZ,
                                   CHAM_desc_t *descZout,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_sgemv_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descZ,
                                   CHAM_desc_t *descZout,
                                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void EXAGEOSTAT_TASK_sdconv(const RUNTIME_option_t *options,
                           int m, int n, int nb,
                           const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAM_desc_t *B, int Bm, int Bn, int ldb);

void EXAGEOSTAT_TASK_dsconv(const RUNTIME_option_t *options,
                           int m, int n, int nb,
                           const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAM_desc_t *B, int Bm, int Bn, int ldb);

void EXAGEOSTAT_MLE_sdmat_reg_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA,
                                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dmse_bivariate_Tile_Async(CHAM_desc_t *descZpre, CHAM_desc_t *descZmiss,
                                            CHAM_desc_t *descserror1, CHAM_desc_t *descserror2, CHAM_desc_t *descserror,
                                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void EXAGEOSTAT_TASK_sexageostat_gemm(const RUNTIME_option_t *options,
                                     CHAM_enum transA, int transB,
                                     int m, int n, int k, int nb,
                                     float alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                                     const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                                     float beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc);

void EXAGEOSTAT_TASK_sexageostat_trsm(const RUNTIME_option_t *options,
                                     CHAM_enum side, CHAM_enum uplo, CHAM_enum transA, CHAM_enum diag,
                                     int m, int n, int nb,
                                     float alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                                     const CHAM_desc_t *B, int Bm, int Bn, int ldb);


void EXAGEOSTAT_MLE_sdmat_register_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA,
                                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_g_to_ng_Tile_Async(CHAM_desc_t *descZ, double* theta,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_ng_loglike_Tile_Async(CHAM_desc_t *descZ, CHAM_desc_t *descsum, double* theta,
                                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_ng_transform_Tile_Async(CHAM_desc_t *descZ, CHAM_desc_t *descflag, const double* theta,
                                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_ng_corr_vec_gen_Tile_Async(CHAM_desc_t *descr, location *l, location *lobs,
                                         double* theta, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dtrace_Tile_Async(CHAM_desc_t *descA, RUNTIME_sequence_t *sequence,
                                    RUNTIME_request_t *request, CHAM_desc_t *descsum, CHAM_desc_t *desctrace);

int EXAGEOSTAT_stride_vec_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descB, CHAM_desc_t *descC,
                                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_tristride_vec_Tile_Async(CHAM_desc_t *descA, CHAM_desc_t *descB, CHAM_desc_t *descC, CHAM_desc_t *descD,
                                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dmloe_mmom_Tile_Async(CHAM_desc_t *descexpr2, RUNTIME_sequence_t *descexpr3,
                                        RUNTIME_sequence_t *descexpr4, RUNTIME_sequence_t *descmloe,
                                        RUNTIME_sequence_t *descmmom, RUNTIME_sequence_t *sequence,
                                        RUNTIME_request_t *request);

int EXAGEOSTAT_MLE_dmse_trivariate_Tile_Async(CHAM_desc_t *descZpre, CHAM_desc_t *descZmiss, CHAM_desc_t *descserror1,
                                             CHAM_desc_t *descserror2, CHAM_desc_t *descserror3,
                                             CHAM_desc_t *descserror, RUNTIME_sequence_t *sequence,
                                             RUNTIME_request_t *request);

#endif /* _EXAGEOSTATCODELETS_H_ */
