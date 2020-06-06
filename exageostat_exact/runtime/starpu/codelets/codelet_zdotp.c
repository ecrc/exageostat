/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_ddotp.c
 *
 * StarPU codelet to Calculate the dot product of the Z vector.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/starpu_exageostat.h"
static void CORE_ddotp_starpu(void *buffers[], void *cl_arg){
    int m, m0;
    int transA, transB;
    int indexC;
    double *A;
    double *B;
    double *C;

    A               = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
    B               = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
    C               = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &m0, &transA, &transB, &indexC);
    if(transA == MorseTrans)
        LAPACKE_sge_trans(LAPACK_COL_MAJOR, m, 1, A, 1, A, 1); 
    if((transB == MorseTrans))
        LAPACKE_sge_trans(LAPACK_COL_MAJOR, m, 1, B, 1, B, 1);
    double local_dot=cblas_ddot(m, A, 1, B, 1);
    C[indexC] += local_dot;
}


static struct starpu_codelet cl_ddotp =
{
    .where          = STARPU_CPU,
    .cpu_funcs      = {CORE_ddotp_starpu},
    .nbuffers       = 3,
    .modes          = {STARPU_R, STARPU_R, STARPU_RW},
    .name           = "ddotp"
};


/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_ddotp_Async- codelet to compute dot product of A.A.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *           Morse descriptor
 *
 * @param[out] descproduct
 *           dot product descriptor.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[in] request
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
int MORSE_MLE_ddotp_Async(MORSE_enum transA, MORSE_enum transB,
        MORSE_desc_t *descA, MORSE_desc_t *descB,
        MORSE_desc_t *descC, MORSE_enum indexC,
        MORSE_sequence_t *sequence, MORSE_request_t *request)

{

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m, m0;
    int tempmm;
    MORSE_desc_t A = *descA;
    MORSE_desc_t B = *descB;
    MORSE_desc_t C = *descC;
    struct starpu_codelet *cl=&cl_ddotp;

    printf("%d - %d -%d\n",descA->m,descB->m,descC->m);
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;
        starpu_insert_task(starpu_mpi_codelet(cl),
                STARPU_VALUE, &tempmm, sizeof(int),
                STARPU_VALUE, &m0,   sizeof(int),
                STARPU_VALUE, &transA,   sizeof(MORSE_enum),
                STARPU_VALUE, &transB,   sizeof(MORSE_enum),
                STARPU_R,  EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, descA->m, descA->n),
                STARPU_R,  EXAGEOSTAT_RTBLKADDR(descB, MorseRealDouble, descB->m, descB->n),
                STARPU_RW, EXAGEOSTAT_RTBLKADDR(descC, MorseRealDouble, descC->m, descC->n),
                STARPU_VALUE, &indexC,   sizeof(MORSE_enum),
                0);

    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    //MORSE_Sequence_Wait(sequence);
    return MORSE_SUCCESS;
}
