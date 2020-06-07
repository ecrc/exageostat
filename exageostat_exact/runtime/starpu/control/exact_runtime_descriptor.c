/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file  exageostat_runtime_descriptor.c
 *
 * StarPU codelets functions header file.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2019-02-13
 *
 **/
#include <stdlib.h>
#include <unistd.h>
#include "../include/starpu_exageostat.h"
static int tag_sep   = 26;
void *EXAGEOSTAT_data_getaddr( const MORSE_desc_t *A, int type, int m, int n )
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t)A->lmt) * nn + mm;

    if (*ptrtile == NULL) {
        int home_node = -1;
        void *user_ptr = NULL;
        int myrank = A->myrank;
        int owner  = A->get_rankof( A, m, n );
        int64_t eltsze = MORSE_Element_Size(type);
        int tempmm = (mm == A->lmt-1) ? (A->lm - mm * A->mb) : A->mb;
        int tempnn = (nn == A->lnt-1) ? (A->ln - nn * A->nb) : A->nb;


        if ( myrank == owner ) {
            user_ptr = A->get_blkaddr(A, m, n);
            if ( user_ptr != NULL ) {
                home_node = STARPU_MAIN_RAM;
            }
        }

        starpu_matrix_data_register( ptrtile, home_node, (uintptr_t) user_ptr,
                BLKLDD(A, m),
                tempmm, tempnn, eltsze );

#ifdef HAVE_STARPU_DATA_SET_COORDINATES
        starpu_data_set_coordinates( *ptrtile, 2, m, n );
#endif

#if defined(CHAMELEON_USE_MPI)
        {
            int64_t block_ind = A->lmt * nn + mm;
            starpu_mpi_data_register(*ptrtile, (A->id << tag_sep) | (block_ind), owner);
        }
#endif /* defined(MORSE_USE_MPI) */
    }

    return *ptrtile;
}
