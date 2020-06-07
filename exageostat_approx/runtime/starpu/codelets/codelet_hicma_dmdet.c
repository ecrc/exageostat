/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dmdet.c
 *
 * StarPU codelet to Calculate determinant of a given triangular matrix (A)
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/starpu_exageostat.h"
static void CORE_dmdet_starpu(void *buffers[],void *cl_arg){
    int m;
    int n;
    double *A;
    int m0;
    int n0;
    double det = 0;
    double *determinant=&det;

    *determinant	= 0;
    A		= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
    determinant	= (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n,  &m0, &n0);

    double local_det=core_dmdet(A, m, n, m0, n0);

    *determinant	+= local_det;
}

static struct starpu_codelet cl_dmdet =
{
    .where		= STARPU_CPU,
    .cpu_funcs	= {CORE_dmdet_starpu},
    .nbuffers	= 2,
    .modes		= {STARPU_R,STARPU_RW},
    .name		= "dmdet"
};

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_dmdet_Tile_Async  - Calculate determinant for triangular matrix.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           descA:  Morse descriptor
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 * @param[out] descdet
 *           descerror:  determinant value
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int HICMA_MLE_dmdet_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request, MORSE_desc_t * descdet) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;

    RUNTIME_options_init(&options, morse, sequence, request);

    int m, m0, n0;
    int tempmm;
    MORSE_desc_t A = *descA;
    struct starpu_codelet *cl=&cl_dmdet;


    for(m=0; m < A.mt; m++)
    {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                STARPU_VALUE, &tempmm,  sizeof(int),
                STARPU_VALUE, &tempmm, sizeof(int),
                STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, m, 0),
                STARPU_VALUE, &m0,   sizeof(int),
                STARPU_VALUE, &n0,   sizeof(int),
                STARPU_RW, EXAGEOSTAT_RTBLKADDR(descdet, MorseRealDouble, 0, 0),
                0);
    }

    //MORSE_TASK_flush_desc( &options, MorseUpperLower, descA);
    //MORSE_TASK_flush_desc( &options, MorseUpperLower, descdet);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    //MORSE_TASK_flush_all();
    //MORSE_TASK_dataflush_all();
    return MORSE_SUCCESS;
}


