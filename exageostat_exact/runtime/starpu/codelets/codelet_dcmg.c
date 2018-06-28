/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file codelet_dcmg.c
 *
 * StarPU codelet to Generate covariance matrix of a set of locations in 2D using Matern kernel.
 *
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include "../include/starpu_exageostat.h"

static void CORE_dcmg_starpu(void *buffers[],void *cl_arg){
        int m, n, m0, n0;
        location * l1;
        location * l2;
        double * theta;
        double * A;
	int distance_metric;
        theta	= (double *) malloc(3* sizeof(double));
        A	= (double *)STARPU_MATRIX_GET_PTR(buffers[0]);

        starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &l1, &l2, &theta[0], &theta[1], &theta[2], &distance_metric);


        core_dcmg(A, m, n, m0, n0, l1, l2, theta, distance_metric);

}

static struct starpu_codelet cl_dcmg =
{
                .where		= STARPU_CPU,
                .cpu_funcs	= {CORE_dcmg_starpu},
                .nbuffers 	= 1,
                .modes		= {STARPU_W},
		.name		= "dcmg"
};

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_MLE_cmg_Tile_Async - Calculate covariance matrix descA.
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
int MORSE_MLE_cmg_Tile_Async(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request, location *l1, location *l2, double * theta , char * dm) {

        MORSE_context_t *morse;
        MORSE_option_t options;
        morse = morse_context_self();

      
  if (sequence->status != MORSE_SUCCESS)
                return -2;
        RUNTIME_options_init(&options, morse, sequence, request);


        int m, n, m0, n0;
	int distance_metric = strcmp(dm,"gc")==0? 1 : 0 ; 
        int tempmm, tempnn;
        MORSE_desc_t A = *descA;
        struct starpu_codelet *cl=&cl_dcmg;


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
                                        STARPU_W, RTBLKADDR(descA, sizeof(double)*ldam*tempnn, m, n),
                                        STARPU_VALUE, &l1, sizeof(location*),
                                        STARPU_VALUE, &l2, sizeof(location*),
                                        STARPU_VALUE, &theta[0], sizeof(double),
                                        STARPU_VALUE, &theta[1], sizeof(double),
                                        STARPU_VALUE, &theta[2], sizeof(double),
                                        STARPU_VALUE, &distance_metric, sizeof(int),
				         0);

                }

        }

        RUNTIME_options_ws_free(&options);
        RUNTIME_options_finalize(&options, morse);
	//MORSE_TASK_dataflush_all();
        //MORSE_TASK_dataflush_all(); is replaced in the new chameleon by  MORSE_Desc_Flush( DESC, sequence );
        MORSE_Desc_Flush( descA, sequence );
 
	return MORSE_SUCCESS;
}



