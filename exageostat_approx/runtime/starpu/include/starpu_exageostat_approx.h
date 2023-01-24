/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file starpu_exageostat_approx.h
 *
 * StarPU codelets functions header file.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#ifndef _STARPU_EXAGEOSTAT_APPROX_H_
#define _STARPU_EXAGEOSTAT_APPROX_H_

#include "../../../cpu_core/include/exageostatcore.h"
#include "chameleon.h"
#include "../../../../include/chameleon_starpu.h"
#include "../../../../include/context.h"
#include "../../../../include/descriptor.h"
#include "../../../../include/chameleon_starpu.h"
#include "starpu_exageostat.h"
#if defined(EXAGEOSTAT_USE_HICMA)
	#include <hicma_struct.h>
#endif
#if defined(CHAMELEON_USE_MPI)
#undef STARPU_REDUX

#define starpu_insert_task starpu_mpi_insert_task

#define starpu_mpi_codelet(_codelet_) MPI_COMM_WORLD, _codelet_
#else

#define starpu_mpi_codelet(_codelet_) _codelet_
#endif

#define EXAGEOSTAT_RTBLKADDR(desc, type, m, n) ( (starpu_data_handle_t)EXAGEOSTAT_data_getaddr( desc, type, m, n ) )

#define mBLKLDD(A, k) A->get_blkldd( A,k )


void *EXAGEOSTAT_data_getaddr(const CHAM_desc_t *A, int type, int m, int n);

int CHAM_MLE_dcmg_diag_Tile_Async(CHAM_enum uplo, CHAM_desc_t *descA, location *l1,
                                  location *l2, location *lm, double* theta,
                                  char* dm, char* kernel_fun, int diag_thick,
                                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void CHAM_TASK_dpotrf_diag(const RUNTIME_option_t *options,
                           CHAM_enum uplo, int n, int nb,
                           const CHAM_desc_t *A, int Am, int An, int lda,
                           int iinfo);
#if defined(EXAGEOSTAT_USE_HICMA)
int EXAGEOSTAT_TLR_MLE_dmdet_Tile_Async(HICMA_desc_t *descA, HICMA_sequence_t *sequence, HICMA_request_t *request,
                               HICMA_desc_t *descdet);

int EXAGEOSTAT_ng_transform_lr_Tile_Async(HICMA_desc_t *descZ, HICMA_desc_t *descflag, const double* theta,
                                      HICMA_sequence_t *sequence, HICMA_request_t *request);

int EXAGEOSTAT_ng_loglike_lr_Tile_Async(HICMA_desc_t *descZ, HICMA_desc_t *descsum, double* theta,
                                    HICMA_sequence_t *sequence, HICMA_request_t *request);
#endif
#endif /* _EXAGEOSTATCODELETS_H_ */
