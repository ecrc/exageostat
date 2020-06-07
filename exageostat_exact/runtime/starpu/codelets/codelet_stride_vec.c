/*
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dzcpy.c
 *
 * StarPU codelet to Copy contents of descriptor to vector
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-03-11
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_stride_vecstarpu(void *buffers[], void *cl_arg){
	int m;
	int tempmm;
	double *A;
        double *B;
        double *C;
	int m0;
	int i=0;
	int j=0;

	A       = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
	B       = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
	C       = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);	
	starpu_codelet_unpack_args(cl_arg,  &tempmm, &m0, &m);
	if(m%2 == 0)
		j=0;
	else
		j=tempmm/2;
	for(i=0;i<tempmm-1;i+=2)
	{
		//printf("%d -%d\n",i, j);
		B[j]=A[i];
		C[j]=A[i+1];
		j++;
	}

}

static struct starpu_codelet cl_stride_vec=
{
	.where          = STARPU_CPU,
	.cpu_funcs      = {CORE_stride_vecstarpu},
	.nbuffers       = 3,
	.modes          = STARPU_R, STARPU_W, STARPU_W,
	.name           = "stride_vec"
};






/***************************************************************************//**
 *
 * @ingroup MORSE_Complex32_t_Tile (single precision).
 *
 *  MORSE_MLE_szcpy_Tile_Async - copy Morse descriptor to vector float*.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[out] descA
 *           Morse descriptor
 *
 * @param[in] r
 *           double*
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
int MORSE_stride_vec_Tile_Async(MORSE_desc_t *descA, MORSE_desc_t *descB, MORSE_desc_t *descC, MORSE_sequence_t *sequence, MORSE_request_t  *request) {
	MORSE_context_t *morse;
	MORSE_option_t options;
	morse = morse_context_self();
	if (sequence->status != MORSE_SUCCESS)
		return -2;
	RUNTIME_options_init(&options, morse, sequence, request);

	int m, m0;
	int tempmm;
	MORSE_desc_t A = *descA;
	MORSE_desc_t B = *descA;
	MORSE_desc_t C = *descA;
	struct starpu_codelet *cl=&cl_stride_vec;
	for (m = 0; m < A.mt; m++) {
		tempmm = m == A.mt-1 ? A.m - m * A.mb : A.mb;
		m0 = m * A.mb;
		starpu_insert_task(starpu_mpi_codelet(cl),
				STARPU_VALUE, &tempmm, sizeof(int),
				STARPU_VALUE, &m0, sizeof(int),
                                STARPU_VALUE, &m, sizeof(int),
				STARPU_R, EXAGEOSTAT_RTBLKADDR(descA, MorseRealDouble, m, 0),
				STARPU_W, EXAGEOSTAT_RTBLKADDR(descB, MorseRealDouble, (int)floor(m/2.0), 0),
				STARPU_W, EXAGEOSTAT_RTBLKADDR(descC, MorseRealDouble, (int)floor(m/2.0), 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
				STARPU_NAME, "stride_vec",
#endif
				0);

	}

	RUNTIME_options_ws_free(&options);
	//MORSE_TASK_dataflush_all();
	return MORSE_SUCCESS;
}
