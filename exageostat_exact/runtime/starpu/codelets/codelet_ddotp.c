/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include "../include/starpu_exageostat.h"
static void CORE_ddotp_starpu(void *buffers[], void *cl_arg){
        int m, m0;
        double * A;
        double * dotproduct;

        dotproduct	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
        A		= (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
        starpu_codelet_unpack_args(cl_arg, &m, &m0);
        double local_dot=cblas_ddot(m, A, 1, A, 1);
        *dotproduct += local_dot;
}

static struct starpu_codelet cl_ddotp =
{
                .where		= STARPU_CPU,
                .cpu_funcs	= {CORE_ddotp_starpu},
                .nbuffers	= 2,
                .modes		= {STARPU_RW,STARPU_R},
		.name		= "ddotp"
};

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_ddotp_Async- colete to compute dot product of A.A.
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
int MORSE_MLE_ddotp_Async(MORSE_desc_t *descA, MORSE_desc_t *descproduct, MORSE_sequence_t *sequence, MORSE_request_t  *request) {

        MORSE_context_t *morse;
        MORSE_option_t options;
        morse = morse_context_self();
        if (sequence->status != MORSE_SUCCESS)
                return -2;
        RUNTIME_options_init(&options, morse, sequence, request);

        int m, m0;
        int tempmm;
        MORSE_desc_t A = *descA;

        struct starpu_codelet *cl=&cl_ddotp;

        for (m = 0; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m - m * A.mb : A.mb;

                m0 = m * A.mb;


                starpu_insert_task(starpu_mpi_codelet(cl),
                                STARPU_VALUE, &tempmm, sizeof(int),
                                STARPU_VALUE, &m0,   sizeof(int),
                                STARPU_RW, RTBLKADDR(descproduct, sizeof(double), 0, 0),
                                STARPU_R, RTBLKADDR(descA, sizeof(double)*tempmm, m, 0),
			         0);
        }
        RUNTIME_options_ws_free(&options);
        RUNTIME_options_finalize(&options, morse);
	//MORSE_TASK_dataflush_all();
        //MORSE_TASK_dataflush_all(); is replaced in the new chameleon by  MORSE_Desc_Flush( DESC, sequence );
        MORSE_Desc_Flush( descA, sequence );
        MORSE_Desc_Flush( descproduct, sequence );
        MORSE_Sequence_Wait(sequence);
        return MORSE_SUCCESS;
}


