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
 * @date 2018-11-11
 *
 **/
#include "../include/starpu_exageostat.h"
static void CORE_dmse_starpu(void *buffers[], void *cl_arg){
    int m, m0, i;
    double * zpre;
    double * zmiss;
    double * serror;
    double local_serror = 0.0;

    serror	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
    zpre	= (double *)STARPU_MATRIX_GET_PTR(buffers[1]);
    zmiss	= (double *)STARPU_MATRIX_GET_PTR(buffers[2]);

    starpu_codelet_unpack_args(cl_arg, &m, &m0);

    for(i = 0; i < m; i++)
    {
        // 		printf("%f, %f, \n",zpre[i], zmiss[i]);
        local_serror += pow((zpre[i]-zmiss[i]), 2);
    }

    *serror += local_serror;
}

static struct starpu_codelet cl_dmse =
{
    .where		= STARPU_CPU,
    .cpu_funcs	= {CORE_dmse_starpu},
    .nbuffers 	= 3,
    .modes		= {STARPU_RW,STARPU_R,STARPU_R},
    .name		= "dmse"
};


//******************************************************************************

static void CORE_smse_starpu(void *buffers[], void *cl_arg){
    int m, m0, i;
    float *zpre;
    float *zmiss;
    float *serror;
    float local_serror = 0.0;

    serror  = (float *)STARPU_MATRIX_GET_PTR(buffers[0]);
    zpre    = (float *)STARPU_MATRIX_GET_PTR(buffers[1]);
    zmiss   = (float *)STARPU_MATRIX_GET_PTR(buffers[2]);

    starpu_codelet_unpack_args(cl_arg, &m, &m0);

    for(i = 0; i < m; i++)
    {
        //              printf("%f, %f, \n",zpre[i], zmiss[i]);
        local_serror += pow((zpre[i]-zmiss[i]), 2);
    }

    *serror += local_serror;
}

static struct starpu_codelet cl_smse =
{
    .where          = STARPU_CPU,
    .cpu_funcs      = {CORE_smse_starpu},
    .nbuffers       = 3,
    .modes          = {STARPU_RW,STARPU_R,STARPU_R},
    .name           = "smse"
};


/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_dmse_Tile_Async - Calculate mean square eror (MSE) scalar value the prediction.
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

int MORSE_MLE_dmse_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss, MORSE_desc_t *descserror, MORSE_sequence_t *sequence, MORSE_request_t  *request){
    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m, m0;
    int tempmm;
    MORSE_desc_t Zpre = *descZpre;
    struct starpu_codelet *cl=&cl_dmse;


    for (m = 0; m < Zpre.mt; m++) {
        tempmm = m == Zpre.mt-1 ? Zpre.m - m * Zpre.mb : Zpre.mb;

        m0 = m * Zpre.mb;

        starpu_insert_task(starpu_mpi_codelet(cl),
                STARPU_VALUE, &tempmm, sizeof(int),
                STARPU_VALUE, &m0,   sizeof(int),
                STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror, MorseRealDouble, 0, 0),
                STARPU_R, EXAGEOSTAT_RTBLKADDR(descZpre, MorseRealDouble, m, 0),
                STARPU_R, EXAGEOSTAT_RTBLKADDR(descZmiss, MorseRealDouble, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dmse",
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




/*******************************************************************************
 *
 * @ingroup MORSE_Complex32_t_Tile  (single precision)
 *
 *  MORSE_MLE_smse_Tile_Async - Calculate mean square eror (MSE) scalar value the prediction.
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

int MORSE_MLE_smse_Tile_Async(MORSE_desc_t *descZpre, MORSE_desc_t *descZmiss, MORSE_desc_t *descserror, MORSE_sequence_t *sequence, MORSE_request_t *request){
    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m, m0;
    int tempmm;
    MORSE_desc_t Zpre = *descZpre;
    struct starpu_codelet *cl=&cl_smse;


    for (m = 0; m < Zpre.mt; m++) {
        tempmm = m == Zpre.mt-1 ? Zpre.m - m * Zpre.mb : Zpre.mb;

        m0 = m * Zpre.mb;

        starpu_insert_task(starpu_mpi_codelet(cl),
                STARPU_VALUE, &tempmm, sizeof(int),
                STARPU_VALUE, &m0,   sizeof(int),
                STARPU_RW, EXAGEOSTAT_RTBLKADDR(descserror, MorseRealFloat, 0, 0),
                STARPU_R, EXAGEOSTAT_RTBLKADDR(descZpre, MorseRealFloat, m, 0),
                STARPU_R, EXAGEOSTAT_RTBLKADDR(descZmiss, MorseRealFloat, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "smse",
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
