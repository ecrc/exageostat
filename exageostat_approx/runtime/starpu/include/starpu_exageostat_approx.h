/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _STARPU_EXAGEOSTAT_APPROX_H_
#define _STARPU_EXAGEOSTAT_APPROX_H_
#include "../../../core/include/exageostatcore.h"
#include "morse.h"
#include "../../../../include/morse_starpu.h"
#include "../../../../include/context.h"
#include "../../../../include/descriptor.h"
#include "../../../../include/chameleon_starpu.h"
#include "starpu_exageostat.h"
//#include "chameleon_starpu.h"
#if defined(CHAMELEON_USE_MPI)
#undef STARPU_REDUX

#define starpu_insert_task starpu_mpi_insert_task

#define starpu_mpi_codelet(_codelet_) MPI_COMM_WORLD, _codelet_
#else

#define starpu_mpi_codelet(_codelet_) _codelet_
#endif

#define EXAGEOSTAT_RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)EXAGEOSTAT_data_getaddr( desc, type, m, n ) )

#define mBLKLDD(A, k) A->get_blkldd( A,k )

//#define A(m,n) (double *)plasma_getaddr(A, m, n)


void *EXAGEOSTAT_data_getaddr( const MORSE_desc_t *A, int type, int m, int n );
int MORSE_MLE_dcmg_diag_Tile_Async(MORSE_enum uplo, MORSE_desc_t *descA, location *l1,
        location *l2, location *lm,  double *theta,
        char *dm, char *kernel_fun, int diag_thick,
        MORSE_sequence_t *sequence, MORSE_request_t  *request);
void MORSE_TASK_dpotrf_diag(const MORSE_option_t *options,
        MORSE_enum uplo, int n, int nb,
        const MORSE_desc_t *A, int Am, int An, int lda,
        int iinfo);
int HICMA_MLE_dmdet_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request, MORSE_desc_t * descdet);

#endif /* _EXAGEOSTATCODELETS_H_ */
