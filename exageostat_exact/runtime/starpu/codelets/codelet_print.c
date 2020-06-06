/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dprint.c
 *
 * StarPU codelet to Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_dprint_starpu(void *buffers[],void *cl_arg){
    int m, n, m0, n0;
    double *A;
    A	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);


    core_dprint(A, m, n, m0, n0);

}

static struct starpu_codelet cl_dprint =
{
    .where		= STARPU_CPU,
    .cpu_funcs	= {CORE_dprint_starpu},
    .nbuffers 	= 1,
    .modes		= {STARPU_RW},
    .name		= "dprint"
};


//******************************************************************************
static void CORE_sprint_starpu(void *buffers[],void *cl_arg){
    int m, n, m0, n0;
    float *A;
    A       = (float *) STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0);

    core_sprint(A, m, n, m0, n0);

}

static struct starpu_codelet cl_sprint =
{
    .where          = STARPU_CPU,
    .cpu_funcs      = {CORE_sprint_starpu},
    .nbuffers       = 1,
    .modes          = {STARPU_RW},
    .name           = "sprint"
};



/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_dprint_Tile_Async - Codelet to generate covariance matrix in descriptor descA in  dense format between two sets of locations (l1, l2) (Matern Kernel).
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *		Upper or lower fill of the matrix.	    
 * 
 * @param[out] descA 
 *		descA:  Morse descriptor that handles the generated covariance matrix.
 *
 * @param[in] sequence
 *		Identifies the sequence of function calls that this call belongs to
 *		(for completion checks and exception handling purposes).
 *
 * @param[out] request
 *		Identifies this function call (for exception handling purposes).
 *
 * @param[in] l1
 *		Location struct of the first input.
 *
 * @param[in] l2
 *		Location struct of the second input.
 *
 * @param[in] theta
 *		Parameter vector that should be used to generate the output covariance matrix.
 *
 * @param[in] dm
 *		Distance metric "euclidean Distance ("ED" -->0) or "Great Circle Distance (GCD) -->1".
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/
int MORSE_MLE_dprint_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();


    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);


    int m, n, m0, n0;
    int tempmm, tempnn;
    MORSE_desc_t A = *descA;
    struct starpu_codelet *cl=&cl_dprint;

    for(m = 0; m < A.mt; m++)
    {
        tempmm = m == A.mt -1 ? A.m- m* A.mb : A.mb;

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt -1 ? A.n - n * A.nb : A.nb;

            m0= m * A.mb;
            n0= n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                    STARPU_VALUE, &tempmm, sizeof(int),
                    STARPU_VALUE, &tempnn, sizeof(int),
                    STARPU_VALUE, &m0, sizeof(int),
                    STARPU_VALUE, &n0, sizeof(int),
                    STARPU_RW, RTBLKADDR(descA, sizeof(double)*ldam*tempnn, m, n),
                    0);

        }

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    return MORSE_SUCCESS;
}


/*******************************************************************************
 *
 * @ingroup MORSE_Complex32_t_Tile
 *
 *  MORSE_MLE_sprint_Tile_Async - Calculate covariance matrix descA.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[out] descA
 *           descA:  Morse descriptor that handles the generated covariance matrix
 *
 * @param[in] l1
 *           location struct of the first input
 *
 * @param[in] l2
 *          location struct of the second input
 *
 * @param[in] theta
 *           parameter vector that should be used to generate the output covariance matrix
 *
 * @param[in] dm
 *           distance metric "euclidean Distance (ED"" or "Great Circle Distance (GCD)"
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 *
 ******************************************************************************/


int MORSE_MLE_sprint_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();


    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);


    int m, n, m0, n0;
    int tempmm, tempnn;
    MORSE_desc_t A = *descA;
    struct starpu_codelet *cl=&cl_sprint;


    for(m = 0; m < A.mt; m++)
    {
        tempmm = m == A.mt -1 ? A.m- m* A.mb : A.mb;

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt -1 ? A.n - n * A.nb : A.nb;

            m0= m * A.mb;
            n0= n * A.nb;
            starpu_insert_task(starpu_mpi_codelet(cl),
                    STARPU_VALUE, &tempmm, sizeof(int),
                    STARPU_VALUE, &tempnn, sizeof(int),
                    STARPU_VALUE, &m0, sizeof(int),
                    STARPU_VALUE, &n0, sizeof(int),
                    STARPU_RW, RTBLKADDR(descA, sizeof(float)*ldam*tempnn, m, n),
                    0);

        }

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    return MORSE_SUCCESS;
}


/*******************************************************************************/
int MORSE_MLE_sdprint_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request, int diag_thick) {

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();


    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);


    int m, n, m0, n0;
    int tempmm, tempnn;
    MORSE_desc_t A = *descA;
    struct starpu_codelet *dcl = &cl_dprint;
    struct starpu_codelet *scl = &cl_sprint;


    for(m = 0; m < A.mt; m++)
    {
        tempmm = m == A.mt -1 ? A.m- m* A.mb : A.mb;

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt -1 ? A.n - n * A.nb : A.nb;

            m0= m * A.mb;
            n0= n * A.nb;

            if(abs(m-n) < diag_thick)
            {
                starpu_insert_task(starpu_mpi_codelet(dcl),
                        STARPU_VALUE, &tempmm, sizeof(int),
                        STARPU_VALUE, &tempnn, sizeof(int),
                        STARPU_VALUE, &m0, sizeof(int),
                        STARPU_VALUE, &n0, sizeof(int),
                        STARPU_RW, RTBLKADDR(descA, sizeof(double)*ldam*tempnn, m, n),
                        0);
            }
            else  //for now is the same because the descriptor should be with one type.
            {

                starpu_insert_task(starpu_mpi_codelet(scl),
                        STARPU_VALUE, &tempmm, sizeof(int),
                        STARPU_VALUE, &tempnn, sizeof(int),
                        STARPU_VALUE, &m0, sizeof(int),
                        STARPU_VALUE, &n0, sizeof(int),
                        STARPU_RW, RTBLKADDR(descA, sizeof(float)*ldam*tempnn, m, n),
                        0);
            }
        }

    }


    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    return MORSE_SUCCESS;
}

