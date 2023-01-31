/**
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
 *                          All rights reserved. 
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014, 2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_dsconv.c
 *
 *  CHAM codelets kernel
 *  CHAM is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAM 1.0.0
 * @author Sameh Abdulah
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-20
 *
 **/
#include "chameleon_starpu.h"
#include "runtime_codelet_d.h"
#include "../include/starpu_exageostat.h"
/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 **/

#if !defined(CHAMELEON_SIMULATION)
static void cl_dsconv_cpu_func(void *descr[], void *cl_arg)
{
	int m;
	int n;
	double  *A;
	float  *B;

	A = (double* )STARPU_MATRIX_GET_PTR(descr[0]);
	B = (float *)STARPU_MATRIX_GET_PTR(descr[1]);
	starpu_codelet_unpack_args(cl_arg, &m, &n);

    core_dsconv(A, B, m, n);
    int i=0;
	for (i=0;i <m ;i++)
	{
        if((float)A[i] != B[i])
		{
            printf("something worng\n");
            printf("%f - %f - %d, %d\n", A[i], B[i], i ,m);
            exit(0);
		}
	}
}

#ifdef CHAMELEON_USE_CUDA
static void cl_dsconv_cuda_func(void *descr[], void *cl_arg)
{
}

#endif /* CHAMELEON_USE_CUDA */
#endif /* !defined(CHAMELEON_SIMULATION) */

static struct starpu_codelet cl_dsconv =
{
	.where          = STARPU_CPU,
	.cpu_funcs      = {cl_dsconv_cpu_func},
	.nbuffers       = 2,
	.modes          = {STARPU_R, STARPU_RW},
	.name           = "dsconv"
};


void EXAGEOSTAT_TASK_dsconv(const RUNTIME_option_t *options,
		int m, int n, int nb,
		const CHAM_desc_t *A, int Am, int An, int lda,
		const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
	(void)nb;
	struct starpu_codelet *codelet = &cl_dsconv;
	void (*callback)(void*) =  NULL;
	int sizeA = lda*m;
	int sizeB = ldb*n;
	int execution_rank = B->get_rankof( B, Bm, Bn );
	int rank_changed=0;
	(void)execution_rank;

	/*  force execution on the rank owning the largest data (tile) */
	int threshold;
	char* env = getenv("CHAMELEON_COMM_FACTOR_THRESHOLD");
	if (env != NULL)
		threshold = (unsigned)atoi(env);
	else
		threshold = 10;
	if ( sizeA > threshold*sizeB ){
		execution_rank = A->get_rankof( A, Am, An );
		rank_changed=1;
	}

	CHAMELEON_BEGIN_ACCESS_DECLARATION;
	CHAMELEON_ACCESS_R(A, Am, An);
	CHAMELEON_ACCESS_RW(B, Bm, Bn);
	if (rank_changed)
		CHAMELEON_RANK_CHANGED(execution_rank);
	CHAMELEON_END_ACCESS_DECLARATION;

	starpu_insert_task(
			starpu_mpi_codelet(codelet),
			STARPU_VALUE,     &m,                  sizeof(int),
			STARPU_VALUE,     &n,                  sizeof(int),
			STARPU_R,         RTBLKADDR(A, double, Am, An),
			//STARPU_VALUE,    &lda,                sizeof(int),
			STARPU_RW,        RTBLKADDR(B, float, Bm, Bn),
			//STARPU_VALUE,    &ldb,                sizeof(int),
			STARPU_PRIORITY,  options->priority,
			STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
			STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
			STARPU_NAME, "dsconv",
#endif
			0);
}




/*
 * Codelet definition
 */
CODELETS(dscov, 2, cl_dsconv_cpu_func, cl_dsconv_cuda_func, STARPU_CUDA_ASYNC)
