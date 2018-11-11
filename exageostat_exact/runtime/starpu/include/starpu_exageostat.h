/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
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
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _EXAGEOSTATCODELETS_H_
#define _EXAGEOSTATCODELETS_H_
#include "../../../core/include/exageostatcore.h"
#include "morse.h"
#include "../../../../include/morse_starpu.h"
#include "../../../../include/context.h"

//#include "chameleon_starpu.h"
#if defined(CHAMELEON_USE_MPI)
#undef STARPU_REDUX

#define starpu_insert_task starpu_mpi_insert_task

#define starpu_mpi_codelet(_codelet_) MPI_COMM_WORLD, _codelet_
#else

#define starpu_mpi_codelet(_codelet_) _codelet_
#endif

#define RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )

#define mBLKLDD(A, k) A->get_blkldd( A,k )

//#define A(m,n) (double *)plasma_getaddr(A, m, n)



//int MORSE_MLE_cmg_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request, location *l1, location *l2, double * theta, char * dm);
int MORSE_MLE_dcmg_Tile_Async(MORSE_enum uplo, MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request, location *l1, location *l2, double *theta , char *dm);
int MORSE_MLE_zcpy_Tile_Async(MORSE_desc_t *descA, double *r, MORSE_sequence_t *sequence, MORSE_request_t  *request);
int MORSE_MLE_dmdet_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request, MORSE_desc_t *descdet);
int MORSE_MLE_ddotp_Async(MORSE_desc_t *descA, MORSE_desc_t *descproduct, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_MLE_dmse_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss, MORSE_desc_t *descserror, MORSE_sequence_t *sequence, MORSE_request_t  *request);
int MORSE_MLE_dgemv_Tile_Async(MORSE_desc_t *descA, MORSE_desc_t *descZ, MORSE_desc_t *descZout, MORSE_sequence_t *sequence, MORSE_request_t  *request);
#endif /* _EXAGEOSTATCODELETS_H_ */
