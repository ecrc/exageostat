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
 * StarPU codelet to Calculate Mean Square Error (MSE) between two vectors.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#include "../include/starpu_exageostat.h"
static void CORE_dmse_bivariate_starpu(void *buffers[], void *cl_arg){
        int m, m0, i;
        double * zpre;
        double * zmiss;
        double * serror1;
        double * serror2;
        double * serror;
        double local_serror1 = 0.0;
        double local_serror2 = 0.0;
        double local_serror  = 0.0;

        serror1	 = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
        serror2  = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
        serror   = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);
        zpre     = (double *)STARPU_MATRIX_GET_PTR(buffers[3]);
        zmiss	 = (double *)STARPU_MATRIX_GET_PTR(buffers[4]);

	starpu_codelet_unpack_args(cl_arg, &m, &m0);

	for(i = 0; i < m; i++)
	{
		if(i%2 == 0)
		{
			local_serror1 += pow((zpre[i]-zmiss[i]), 2);
		//	printf("===%f\n", local_serror1);
		}		else
			local_serror2 += pow((zpre[i]-zmiss[i]), 2);
		local_serror += pow((zpre[i]-zmiss[i]), 2);

	//	printf("%f, %f, %f\n", local_serror1, local_serror2, local_serror);
	}

	*serror1 += local_serror1;
	*serror2 += local_serror2;
	*serror  += local_serror;
}

static struct starpu_codelet cl_dmse_bivariate =
{
	.where		= STARPU_CPU,
	.cpu_funcs	= {CORE_dmse_bivariate_starpu},
	.nbuffers 	= 5,
	.modes		= {STARPU_RW, STARPU_RW, STARPU_RW, STARPU_R, STARPU_R},
	.name		= "dmse_bivariate"
};



/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_dmse_bivariate_Tile_Async - Calculate mean square eror (MSE) scalar value the prediction.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] descZpre
 *           descZpre:  Observed measurements descZpre
 *
 * @param[in] descZmiss
 *           descZpre: Missing measurements descZpre
 *
 * @param[out] descserror
 *           descerror:  Mean Square Error (MSE)
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

int MORSE_MLE_dmse_bivariate_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss, MORSE_desc_t *descserror1,  MORSE_desc_t *descserror2,  MORSE_desc_t *descserror,  MORSE_sequence_t *sequence, MORSE_request_t  *request){
	MORSE_context_t *morse;
	MORSE_option_t options;
	morse = morse_context_self();
	if (sequence->status != MORSE_SUCCESS)
		return -2;
	RUNTIME_options_init(&options, morse, sequence, request);

	int m, m0;
	int tempmm;
	MORSE_desc_t Zpre = *descZpre;
	struct starpu_codelet *cl=&cl_dmse_bivariate;


	for (m = 0; m < Zpre.mt; m++) {
		tempmm = m == Zpre.mt-1 ? Zpre.m - m * Zpre.mb : Zpre.mb;

		m0 = m * Zpre.mb;

		starpu_insert_task(starpu_mpi_codelet(cl),
				STARPU_VALUE, &tempmm, sizeof(int),
				STARPU_VALUE, &m0,   sizeof(int),
				STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror1, MorseRealDouble, 0, 0),
				STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror2, MorseRealDouble, 0, 0),
				STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror, MorseRealDouble, 0, 0),
				STARPU_R, EXAGEOSTAT_RTBLKADDR(descZpre, MorseRealDouble, m, 0),
				STARPU_R, EXAGEOSTAT_RTBLKADDR(descZmiss, MorseRealDouble, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
				STARPU_NAME, "dmse_bivariate",
#endif
				0);
	}

	//MORSE_TASK_flush_desc( &options, MorseUpperLower, descZpre);
	//MORSE_TASK_flush_desc( &options, MorseUpperLower, descZmiss);
	//MORSE_TASK_flush_desc( &options, MorseUpperLower, descserror);
	RUNTIME_options_ws_free(&options);
	RUNTIME_options_finalize(&options, morse);
	//MORSE_TASK_dataflush_all();
	MORSE_Sequence_Wait(sequence);
	return MORSE_SUCCESS;
}


