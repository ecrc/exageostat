/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_misc.c
 *
 * Auxiliary functions that are needed by ExaGeoStat.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#include "../include/MLE_misc.h"
//***************************************************************************************

void init_data_values(MLE_data *data)
    //! initiate data struct with default values
    /*!
     * @param[in] data: MLE_data struct with different MLE inputs.
     * */
{
    data->variance 			        = 1.0;		///< Variance parameter.
    data->computation 		        = "exact";	///< exact or approx computation.
    data->c_fun	                    = "matern"; ///< exact or approx computation.
    data->test   		            = 0;      	///< test or real-data computation.
    data->iter_count        	    = 0;  		///< number of iterations to converge.
    data->l1.x                      = NULL;   	///< 2D locations for the first dataset (X vector).
    data->l1.y			            = NULL;		//< 2D locations for the first dataset (Y vector).
    data->l1.z		         	    = NULL;     //< 3D locations for the first dataset (Y vector)
    data->lmiss.x           	    = NULL;		///< 2D locations for the missing data (X vector) (prediction stage).
    data->lmiss.y             	    = NULL;		///< 2D locations for the missing data (Y vector) (prediction stage).
    data->lmiss.z                   = NULL;     ///< 2D locations for the missing data (Y vector) (prediction stage).
    data->lobs.x            	    = NULL;		///< 2D locations for the observed data (X vector) (prediction stage).
    data->lobs.y             	    = NULL;		///< 2D locations for the observed data (Y vector) (prediction stage).
    data->lm.x            		    = NULL; 	///< 2D locations for the median data point (X vector).
    data->lm.y                      = NULL;     ///< 2D locations for the median data point (Y vector).
    data->lobs.z                    = NULL;     ///< 2D locations for the observed data (Y vector) (prediction stage).
    data->descC               	    = NULL;		///< Covariance matrix C descriptor.
    data->descZ             	    = NULL;		///< Measurements Z descriptor.
    data->Adense             	    = NULL;		///< Dense matrix descriptor in the case of approximation mode - accuracy check.
    data->Adense2              	    = NULL;		///< Dense matrix descriptor2 in the case of approximation mode - accuracy check.
    data->descZcpy                  = NULL;		///< A copy of Measurements Z descriptor.
    data->descdet                   = NULL;		///< Determinant descriptor.
    data->descproduct               = NULL;		///< Dot product descriptor.
    data->descZmiss                 = NULL;		///< Missing measurements descriptor.
    data->descC12                   = NULL;		///< Covariance Matrix C12 descriptor.
    data->descC22                   = NULL;		///< Covariance Matrix C22 descriptor.
    data->descZactual               = NULL;		///< Actual Measurements Z descriptor.
    data->descZobs                  = NULL;		///< observed Measurements Z descriptor.
    data->descmse                   = NULL;		///< Mean Square Error (MSE) descriptor.
    data->sequence                  = NULL;		///< MORSE sequence.
    data->request                   = NULL;		///< MORSE request.
    data->verbose                   = 0;		///< Verbose indicator.
    data->log                       = 0;		///< Log files generation indicator, 0-->no, 1-->yes.
    data->avg_exec_time_per_iter    = 0;		///< Avergae execution time per iteration (only used in verbose mode).
    data->total_exec_time           = 0.0;		///< Total execution time (only used in verbose mode).
    data->avg_flops_per_iter        = 0.0;		///< Avergae flops per iteration (only used in verbose mode).
    data->final_loglik              = 0.0;		///< Final log likelihood value.
    data->locsFPath                 = "";		///< Locations file path -- in the case of real dataset (real mode).
    data->obsFPath                  = "";		///< Observations file path --  in the case of real dataset (real mode).
    data->obsFPath2                 = "";       ///< Observations file path --  in the case of real dataset (real mode).
    data->actualZFPath              = "";		///< Actual observations file path -- in the case of prediction.
    data->actualZFPath2             = "";       ///< Actual observations file path -- in the case of prediction.
    data->actualZLocFPath           = "";		///< Actial locations file path -- in the case of prediction.
    data->det                       = 0.0;		///< determinant value.
    data->dotp                      = 0.0;		///< dot product value.
    data->mserror                   = 0.0;		///< Mean Square Error (MSE) value.
    data->mserror1                  = 0.0;      ///< Mean Square Error (MSE) value.
    data->mserror2                  = 0.0;      ///< Mean Square Error (MSE) value.
    data->dm                        = "ed";		///< Distance metric to be used ed->Euclidian Distance -- gcd->Great Circle Distance.
    data->diag_thick                = 0;		///< The thick of used diagonal in the case of diagonal approximation approach.
    data->nFileLog                  = NULL;		///< log file name (only used if log -->1).
    data->pFileLog                  = NULL;		///< log file path (only used if log -->1).
    data->hicma_maxrank             = 0;		///< Max Rank in the case of LR-HiCMA approx
    data->hicma_data_type           = 0;		/// to define the problem typr to HiCMA (HICMA_STARSH_PROB_GEOSTAT (Synthetic) or HICMA_STARSH_PROB_GEOSTAT_POINT (real))
    data->hicma_descC               = NULL;		///< HiCMA descC descriptor (for accuracy check).
    data->hicma_descZ               = NULL;		///< HiCMA descZ descriptor.
    data->hicma_descCD              = NULL;		///< HiCMA descCD descriptor.
    data->hicma_descCUV             = NULL;		///< HiCMA descCUV descriptor.
    data->hicma_descCrk             = NULL;		///< HiCMA descCrk descriptor.
    data->hicma_descZcpy            = NULL;		///< A copy of Measurements Z descriptor.
    data->hicma_descdet             = NULL;		///< Determinant descriptor.
    data->hicma_descproduct         = NULL;		///< Dot product descriptor.
    data->hicma_descC12D        	= NULL;		///< HiCMA descCD descriptor.
    data->hicma_descC12UV        	= NULL;		///< HiCMA descCUV descriptor.
    data->hicma_descC12rk        	= NULL;		///< HiCMA descCrk descriptor.
    data->hicma_descC22D        	= NULL;		///< HiCMA descCD descriptor.
    data->hicma_descC22UV        	= NULL;		///< HiCMA descCUV descriptor.
    data->hicma_descC22rk           = NULL;		///< HiCMA descCrk descriptor.
    data-> hicma_acc                = 9; 		///< Accuracy in the case of LR-HiCMA approx.
    data->hsequence                 = NULL;  	///< HiCMA sequence.
    data->hrequest                  = NULL;   	///< HiCMA request.
    data->opt_tol			        = 5;        ///< The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
    data->opt_max_iters		        = -1;     	///< Maximum number of mle iterations.
    data->ooc                       = 0;        ///< Support Out-Of-Core execution, 0-->no, 1-->yes.
    data->kernel_fun         	    = "univariate_matern_stationary";///< stationary_matern, or non_stationary_matern.
    data->precision                 = 0;	///< Double, single, or mixed.
    data->desck_t 			        = NULL;
    data->desck_a 			        = NULL;
    data->desck_ttmp                = NULL;
    data->desck_atmp                = NULL;
    data->descK_ttmp 		        = NULL;
    data->descK_t 			        = NULL;
    data->descK_a 			        = NULL;
    data->descexpr1 		        = NULL;
    data->descexpr2 		        = NULL;
    data->descexpr3 		        = NULL;
    data->descexpr4 		        = NULL;
    data->expr1 			        = 0;
    data->expr2 			        = 0;
    data->expr3 			        = 0;
    data->expr4 			        = 0;
    data->mloe_mmom 		        = 0;
    data->mloe_mmom_async           = 0;
    data->mloe	                    = 0;
    data->mmom                      = 0;
    data->checkpoint_file		    = "";
    data->recovery_file		        = "";

    //init result struct
    results.problem_size		    = -1;
    results.computation		        = NULL;
    results.kernel		        	= NULL;
    results.ds_type			        = NULL;
    results.precision	        	= NULL;
    results.z_sample	        	= -1;
    results.dense_ts		        = -1;
    results.lr_ts			        = -1;
    results.lr_acc                  = -1;
    results.lr_maxrank              = -1;
    results.ncores			        = -1;
    results.ngpus			        = -1;
    results.p			            = -1;
    results.q			            = -1;
    results.num_params		        = -1;
    results.initial_theta		    = NULL;
    results.starting_theta		    = NULL;
    results.estimated_theta		    = NULL;
    results.final_loglik		    = -1;
    results.time_per_iteration	    = -1;
    results.flops_per_iteration	    = -1;
    results.total_mle_time		    = -1;
    results.mse_pred1		        = -1;
    results.mse_pred2		        = -1;
    results.mse_pred			    = -1;
    results.total_pred_time		    = -1;
    results.total_pred_flops		= -1;
    results.mloe			        = -1;
    results.mmom			        = -1;
    results.mloe_exec		        = NULL;
    results.total_mloe_mmom_time	= -1;
    results.total_mloe_mmom_flops	= -1;
    data->checkpoint_file		    = "";
    data->recovery_file		        = "";
    data->time_slots		        = 0;
    data->idw			            = 0 ;
}




static uint32_t Compact1By1(uint32_t x)
    //! Collect every second bit into lower part of input
{
    x &= 0x55555555;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >>  1)) & 0x33333333;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >>  2)) & 0x0f0f0f0f;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >>  4)) & 0x00ff00ff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >>  8)) & 0x0000ffff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

static uint32_t DecodeMorton2X(uint32_t code)
    //! Decode first input
{
    return Compact1By1(code >> 0);
}

static uint32_t DecodeMorton2Y(uint32_t code)
    //! Decode second input
{
    return Compact1By1(code >> 1);
}

static uint32_t Part1By1(uint32_t x)
    //! Spread lower bits of input
{
    x &= 0x0000ffff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x <<  8)) & 0x00ff00ff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x <<  4)) & 0x0f0f0f0f;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x <<  2)) & 0x33333333;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x <<  1)) & 0x55555555;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}
static int compare_uint32(const void *a, const void *b)
    //! Compare two uint32_t
{
    uint32_t _a = *(uint32_t *)a;
    uint32_t _b = *(uint32_t *)b;
    if(_a < _b) return -1;
    if(_a == _b) return 0;
    return 1;
}


static int compare_sdata(const void *a, const void *b)
    //! Compare two sdata
{
    sdata _a = *(sdata *)a;
    sdata _b = *(sdata *)b;
    if(_a.xy < _b.xy) return -1;
    if(_a.xy == _b.xy) return 0;
    return 1;
}


static uint32_t EncodeMorton2(uint32_t x, uint32_t y)
    //! Encode two inputs into one
{
    return (Part1By1(y) << 1) + Part1By1(x);
}

static void zsort_locations(int n, location * locations)
//! Sort in Morton order (input points must be in [0;1]x[0;1] square])
{
    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y;
    uint32_t z[n];
    // Encode data into vector z
    for(i = 0; i < n; i++)
    {
        x = (uint16_t)(locations->x[i]*(double)UINT16_MAX +.5);
        y = (uint16_t)(locations->y[i]*(double)UINT16_MAX +.5);
        //printf("%f %f -> %u %u\n", points[i], points[i+n], x, y);
        z[i] = EncodeMorton2(x, y);
    }
    // Sort vector z
    qsort(z, n, sizeof(uint32_t), compare_uint32);
    // Decode data from vector z
    for(i = 0; i < n; i++)
    {
        x = DecodeMorton2X(z[i]);
        y = DecodeMorton2Y(z[i]);
        locations->x[i] = (double)x/(double)UINT16_MAX;
        locations->y[i] = (double)y/(double)UINT16_MAX;
        //printf("%lu (%u %u) -> %f %f\n", z[i], x, y, points[i], points[i+n]);
    }
}

static uint64_t Part1By3(uint64_t x)
    // Spread lower bits of input
{
    x &= 0x000000000000ffff;
    // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x << 24)) & 0x000000ff000000ff;
    // x = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    x = (x ^ (x << 12)) & 0x000f000f000f000f;
    // x = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    x = (x ^ (x << 6)) & 0x0303030303030303;
    // x = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    x = (x ^ (x << 3)) & 0x1111111111111111;
    // x = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    return x;
}
static uint64_t EncodeMorton3(uint64_t x, uint64_t y, uint64_t z)
    // Encode 3 inputs into one
{
    return (Part1By3(z) << 2) + (Part1By3(y) << 1) + Part1By3(x);
}

static uint64_t Compact1By3(uint64_t x)
    // Collect every 4-th bit into lower part of input
{
    x &= 0x1111111111111111;
    // x = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    x = (x ^ (x >> 3)) & 0x0303030303030303;
    // x = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    x = (x ^ (x >> 6)) & 0x000f000f000f000f;
    // x = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    x = (x ^ (x >> 12)) & 0x000000ff000000ff;
    // x = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    x = (x ^ (x >> 24)) & 0x000000000000ffff;
    // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}
static uint64_t DecodeMorton3X(uint64_t code)
    // Decode first input
{
    return Compact1By3(code >> 0);
}

static uint64_t DecodeMorton3Y(uint64_t code)
    // Decode second input
{
    return Compact1By3(code >> 1);
}

static uint64_t DecodeMorton3Z(uint64_t code)
    // Decode third input
{
    return Compact1By3(code >> 2);
}

static void zsort_locations_3d(int n, location * locations)
    //! Sort in Morton order (input points must be in [0;1]x[0;1] square])
{
    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y, z;
    uint64_t Z[n];
    // Encode data into vector z
    for(i = 0; i < n; i++)
    {
        x = (uint16_t)(locations->x[i]*(double)UINT16_MAX +.5);
        y = (uint16_t)(locations->y[i]*(double)UINT16_MAX +.5);
        z = (uint16_t)(locations->z[i]*(double)UINT16_MAX +.5);
        Z[i] = EncodeMorton3(x, y, z);
    }
    // Sort vector z
    qsort(Z, n, sizeof(uint64_t), compare_uint32);
    // Decode data from vector z
    for(i = 0; i < n; i++)
    {
        x = DecodeMorton3X(Z[i]);
        y = DecodeMorton3Y(Z[i]);
        z = DecodeMorton3Z(Z[i]);
        locations->x[i] = (double)x/(double)UINT16_MAX;
        locations->y[i] = (double)y/(double)UINT16_MAX;
        locations->z[i] = (double)z/(double)UINT16_MAX;
    }
}


static void radix_sort_recursive(uint32_t *data, int count, int ndim,
        int *order, int *tmp_order, int sdim, int sbit,
        int lo, int hi)
    // Hierarchical radix sort to get Z-order of particles.
    // This function is static not to be visible outside this module.
{
    int i, lo_last = lo, hi_last = hi;
    uint32_t *sdata = data+sdim*count;
    uint32_t check = 1 << sbit;
    for(i = lo; i <= hi; i++)
    {
        if((sdata[order[i]] & check) == 0)
        {
            tmp_order[lo_last] = order[i];
            lo_last++;
        }
        else
        {
            tmp_order[hi_last] = order[i];
            hi_last--;
        }
    }
    for(i = lo; i <= hi; i++)
        order[i] = tmp_order[i];
    if(sdim > 0)
    {
        if(lo_last-lo > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, sdim-1,
                    sbit, lo, lo_last-1);
        if(hi-hi_last > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, sdim-1,
                    sbit, hi_last+1, hi);
    }
    else if(sbit > 0)
    {
        if(lo_last-lo > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, ndim-1,
                    sbit-1, lo, lo_last-1);
        if(hi-hi_last > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, ndim-1,
                    sbit-1, hi_last+1, hi);
    }
}

static int radix_sort(uint32_t *data, int count, int ndim,
        int *order)
    // Auxiliary sorting function for starsh_particles_zsort_inpace().
    // This function is static not to be visible outside this module.
{
    int *tmp_order = (int *) malloc(count * sizeof(int));
    radix_sort_recursive(data, count, ndim, order, tmp_order, ndim-1, 31, 0, count-1);
    free(tmp_order);
    return 0;
}

int locations_obs_zsort_inplace(int n, location *locations, double *z)
    //! Sort particles in Z-order (Morton order).
    /*! This function must be used after initializing @ref STARSH_particles with
     * your own data by starsh_particles_init() or starsh_particles_new().
     *
     * @sa starsh_particles_init(), starsh_particles_new().
     * @ingroup app-particles
     * */
{
    int i;
    int j;//new_j, tmp_j;
    int count = n;
    int ndim = 2;
    int info;
    double *point = (double *) malloc(ndim * count * sizeof(double));
    double *ptr1;
    double *minmax; // min is stored in lower part, max is stored in upper part
    //double tmp_x;
    for(i = 0; i< count; i++)	    
    {
        point[i]       = locations->x[i];
        point[i+count] = locations->y[i];
    }

    // I think count should be replaced by ndim 
    minmax = (double *) malloc(2 * count * sizeof(double));

    for(i = 0; i < ndim; i++)
    {
        ptr1 = point+i*count; // i-th dimension
        minmax[i] = ptr1[0];
        minmax[i+ndim] = minmax[i];
        for(j = 1; j < count; j++)
        {
            if(minmax[i] > ptr1[j])
                minmax[i] = ptr1[j];
            else if (minmax[i+ndim] < ptr1[j])
                minmax[i+ndim] = ptr1[j];
        }
    }
    // Now minmax[0:ndim] and minmax[ndim:2*ndim] store minimal and maximal
    // values of coordinates
    uint32_t *uint_point = (uint32_t *) malloc(ndim * count * sizeof(uint32_t));
    uint32_t *uint_ptr1;
    double min, range;
    for(i = 0; i < ndim; i++)
    {
        uint_ptr1 = uint_point+i*count;
        ptr1 = point+i*count;
        min = minmax[i];
        range = minmax[i+ndim]-min;
        for(j = 0; j < count; j++)
            uint_ptr1[j] = (ptr1[j]-min)/range*UINT32_MAX;
    }
    free(minmax);
    // Now uint_ptr1 contains initial coordinates, rescaled to range
    // [0, UINT32_MAX] and converted to uint32_t type to use special radix sort
    // Prepare indexes to store sort order
    int *order = (int *) malloc(count * sizeof(int));
    for(j = 0; j < count; j++)
        order[j] = j;
    info = radix_sort(uint_point, count, ndim, order);
    if(info != 0)
    {
        free(uint_point);
        free(order);
        return info;
    }
    double *new_point = (double *) malloc(ndim * count * sizeof(double));
    double *new_z     = (double *) malloc( count * sizeof(double));
    for(j = 0; j < count; j++)
    {
        for(i = 0; i < ndim; i++)
            new_point[count*i+j] = point[count*i+order[j]];
        new_z[j] = z[order[j]];
    }
    for(i = 0; i< count; i++)
    {
        locations->x[i] = new_point[i];
        locations->y[i] = new_point[i+count];
        z[i] = new_z[i];
        //printf("%f - %f - %f\n",locations->x[i], locations->y[i], z[i]);

    }
    free(new_point);
    free(point);
    free(new_z);
    free(uint_point);
    free(order);
    return 0;
}



//int locations_obs_zsort_inplace(int n, location *locations, double *z)
//! Sort particles in Z-order (Morton order).
/*! This function must be used after initializing @ref STARSH_particles with
 * your own data by starsh_particles_init() or starsh_particles_new().
 *
 * @sa starsh_particles_init(), starsh_particles_new().
 * @ingroup app-particles
 * */
/*{
  int i;
  int j, new_j, tmp_j;
  int count = n;
  int ndim = 2;
  int info;
  double *point = data->point;  //********
  double *ptr1;
  double *minmax; // min is stored in lower part, max is stored in upper part
  double tmp_x;
  minmax = (double *) malloc(2 * count * sizeof(double));

//Stop
for(i = 0; i < ndim; i++)
{
ptr1 = point+i*count; // i-th dimension
minmax[i] = ptr1[0];
minmax[i+ndim] = minmax[i];
for(j = 1; j < count; j++)
{
if(minmax[i] > ptr1[j])
minmax[i] = ptr1[j];
else if (minmax[i+ndim] < ptr1[j])
minmax[i+ndim] = ptr1[j];
}
}
// Now minmax[0:ndim] and minmax[ndim:2*ndim] store minimal and maximal
// values of coordinates
uint32_t *uint_point;
STARSH_MALLOC(uint_point, count*ndim);
uint32_t *uint_ptr1;
double min, range;
for(i = 0; i < ndim; i++)
{
uint_ptr1 = uint_point+i*count;
ptr1 = point+i*count;
min = minmax[i];
range = minmax[i+ndim]-min;
for(j = 0; j < count; j++)
uint_ptr1[j] = (ptr1[j]-min)/range*UINT32_MAX;
}
free(minmax);
// Now uint_ptr1 contains initial coordinates, rescaled to range
// [0, UINT32_MAX] and converted to uint32_t type to use special radix sort
// Prepare indexes to store sort order
STARSH_int *order;
STARSH_MALLOC(order, count);
for(j = 0; j < count; j++)
order[j] = j;
info = radix_sort(uint_point, count, ndim, order);
if(info != STARSH_SUCCESS)
{
free(uint_point);
free(order);
return info;
}
double *new_point;
STARSH_MALLOC(new_point, count*ndim);
for(j = 0; j < count; j++)
{
for(i = 0; i < ndim; i++)
new_point[count*i+j] = point[count*i+order[j]];
}
data->point = new_point;
free(point);
free(uint_point);
free(order);
return STARSH_SUCCESS;
}
*/

void zsort_locations_obs(int n, location *locations, double *z)
    //! Sort in Morton order (input points must be in [0;1]x[0;1] square])
{
    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y;
    sdata z_struct[n];

    // Encode data into vector z
    for(i = 0; i < n; i++)
    {
        x = (uint16_t)(locations->x[i]*(double)UINT16_MAX +.5);
        y = (uint16_t)(locations->y[i]*(double)UINT16_MAX +.5);
        z_struct[i].xy = EncodeMorton2(x, y);
        z_struct[i].z  = z[i];
    }

    // Sort vector z
    qsort(z_struct, n, sizeof(sdata), compare_sdata);
    // Decode data from vector z
    for(i = 0; i < n; i++)
    {
        x = DecodeMorton2X(z_struct[i].xy);
        y = DecodeMorton2Y(z_struct[i].xy);
        locations->x[i] = (double)x/(double)UINT16_MAX;
        locations->y[i] = (double)y/(double)UINT16_MAX;
        z[i]		= z_struct[i].z;
    }
}


double uniform_distribution(double rangeLow, double rangeHigh) 
    //! Generate uniform distribution between rangeLow , rangeHigh
{
    // unsigned int *seed = &exageostat_seed;
    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = rangeHigh - rangeLow;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

location* GenerateXYLoc(int n, int seed) 
    //! Generate XY location for exact computation (MOORSE)        
{
//initalization
    int i = 0, index = 0, j = 0;
    srand(seed);
    location *locations = (location *) malloc(sizeof(location *));
    //Allocate memory
    locations->x = (double* ) malloc(n * sizeof(double));
    locations->y = (double* ) malloc(n * sizeof(double));
    locations->z = NULL;

    int sqrtn = ceil(sqrt(n));

    int *grid = (int *) calloc((int) sqrtn, sizeof(int));

    for (i = 0; i < sqrtn; i++) {
        grid[i] = i + 1;
    }

    for (i = 0; i < sqrtn && index < n; i++)
        for (j = 0; j < sqrtn && index < n; j++) {
            locations->x[index] = (grid[i] - 0.5 + uniform_distribution(-0.4, 0.4)) / sqrtn;
            locations->y[index] = (grid[j] - 0.5 + uniform_distribution(-0.4, 0.4)) / sqrtn;
            index++;
        }
    free(grid);
    zsort_locations(n, locations);
    return locations;
}


location* GenerateXYZLoc(int n, int seed)
    //! Generate XY location for exact computation (MOORSE)
{
    //initalization
    int i = 0 ,index = 0, j = 0, k = 0;
    // unsigned int *seed = &exageostat_seed;
    srand(seed);
    location* locations = (location *) malloc( sizeof(location*));
    //Allocate memory
    locations->x            = (double *) malloc(n * sizeof(double));
    locations->y            = (double *) malloc(n * sizeof(double));
    locations->z            = (double *) malloc(n * sizeof(double));


    int cbrtn = cbrt(n);

    //Check if the input is square number or not
    if(pow(cbrtn, 3) != n)
    {
        printf("Please use a perfect cubic number to generate a valid synthetic dataset.....\n\n");
        exit(0);
    }

    int *grid = (int *) calloc((int)cbrtn, sizeof(int));
    for(i = 0; i < cbrtn; i++)
    {
        grid[i] = i+1;
    }

    printf("cbrtn:%d\n", cbrtn);

    for(i = 0; i < cbrtn; i++)
        for(j = 0; j < cbrtn; j++)
            for(k = 0; k < cbrtn; k++)
            {
                locations->x[index] = (grid[i]-0.5+uniform_distribution(-0.4, 0.4))/cbrtn;
                locations->y[index] = (grid[j]-0.5+uniform_distribution(-0.4, 0.4))/cbrtn;
                locations->z[index] = (grid[k]-0.5+uniform_distribution(-0.4, 0.4))/cbrtn;
                index++;
            }
    free(grid);
    zsort_locations_3d(n, locations);
    return locations;
}


location* GenerateXYLoc_ST(int n, int t_slots, int seed)
    //! Generate XY location for exact computation (MOORSE)
{
    //initalization
    int i = 0 ,index = 0, j = 0;
    n/=t_slots; ///to be modified
    // unsigned int *seed = &exageostat_seed;
    srand(seed);
    location* locations = (location *) malloc( sizeof(location*));
    //Allocate memory
    locations->x            = (double *) malloc(n * t_slots * sizeof(double));
    locations->y            = (double *) malloc(n * t_slots * sizeof(double));
    locations->z            = (double *) malloc(n * t_slots * sizeof(double));
    // if(strcmp(locs_file, "") == 0)
    // {
    int sqrtn = sqrt(n);

    //Check if the input is square number or not
    if(pow(sqrtn,2) != n)
    {
        printf("n=%d, Please use a perfect square number to generate a valid synthetic dataset.....\n\n", n);
        exit(0);
    }

    int *grid = (int *) calloc((int)sqrtn, sizeof(int));

    for(i = 0; i < sqrtn; i++)
    {
        grid[i] = i+1;
    }

    for(i = 0; i < sqrtn; i++)
        for(j = 0; j < sqrtn; j++){
            locations->x[index] = (grid[i]-0.5+uniform_distribution(-0.4, 0.4))/sqrtn;
            locations->y[index] = (grid[j]-0.5+uniform_distribution(-0.4, 0.4))/sqrtn;
            //printf("%f, %f\n", locations->x[index], locations->y[index]);
            index++;
        }

    for(j=0;j<t_slots;j++)
        for(i=0;i< n;i++)
        {
            locations->x[i+j*n]= locations->x[i]; 
            locations->y[i+j*n]= locations->y[i];
            locations->z[i+j*n] = (double)j+1;
            //				printf("%f, %f, %f\n", locations->x[i+j*n], locations->y[i+j*n], locations->z[i+j*n]);
            //printf("%f, %f, %f\n", locations->x[0], locations->y[0], locations->z[0]);
            //printf("%d\n",i+j*n);
        }

    //	exit(0);
    //	free(grid);
    //	zsort_locations_3d(n, locations);
    //printf("%f, %f, %f\n", locations->x[18], locations->y[18], locations->z[18]);
    //exit(0);
    //printf("%d\n",i+j*n);
    return locations;
}

void print_dmatrix(char* desc, int m, int n, double* a, int lda)
    //! print matrix contents (only for testing accuracy should be removed in release)
{
    int i, j;
    fprintf(stderr,"\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            fprintf(stderr, " %6.4e", a[i + lda * j]);
        fprintf(stderr,"\n");
    }
}


void print_diagonal(char* desc, int m, double* a, int lda)
//! print matrix contents (only for testing accuracy should be removed in release)
{
    int i, j;
    fprintf(stderr,"\n %s\n", desc);
    for (i = 0; i < m; i++) {
        fprintf(stderr, " %6.4e - ", a[i + lda * i]);
    }
    fprintf(stderr,"\n");
}


void print_smatrix(char* desc, int m, int n, float* a, int lda)
//! print matrix contents (only for testing accuracy should be removed in release)
{
    int i, j;
    fprintf(stderr,"\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            fprintf(stderr, " %6.4e", a[i + lda * j]);
        fprintf(stderr,"\n");
    }
}


int countlines(char *filename)
    //!count the number of samples in real running mode
{
    FILE *fp = fopen(filename,"r");
    int ch = 0;
    int lines = 0;

    if (fp == NULL)
    {
        fprintf(stderr,"cannot open locations file\n");
        return 0;
    }

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
            lines++;
    }

    fclose(fp);

    //Excluding header line
    return (lines);
}


void write_to_file(char * path, int matrix_size,int ncores,int tile_size, int test, char *computation, int async, char *obsFPath,double total_exec_time,double avg_exec_time_per_iter, double avg_flops_per_iter , int p_grid, int q_grid, double final_loglik, int n)
//! write results in  detail (only for testing accuracy should be removed in release)
{
    FILE *pFile;
    double peakperfomrance;
    double percent;    
    pFile = fopen(path,"a");
    if(pFile == NULL) {
        fprintf(stderr,"Cannot access the results path(1)\n");
        exit(0);
    }
    peakperfomrance = q_grid*p_grid*ncores*16*2.3;
    percent = avg_flops_per_iter/peakperfomrance;

    fprintf(pFile, "%d\t", n);
    fprintf(pFile, "%d\t", ncores);
    fprintf(pFile, "%d\t", q_grid*p_grid);
    fprintf(pFile, "%f\t", total_exec_time);
    fprintf(pFile, "%f\t", avg_exec_time_per_iter);
    fprintf(pFile, "%f\t\t", avg_flops_per_iter);
    fprintf(pFile, "%f\t", peakperfomrance);
    fprintf(pFile, "%f\t", percent);
    fprintf(pFile, "%d-", tile_size);
    //	fprintf(pFile, "%s-", ikernel);
    fprintf(pFile, "%s-", computation);
    if(async == 0)
        fprintf(pFile, "SYNC-");
    else
        fprintf(pFile, "ASYNC-");
    if(test == 1)
        fprintf(pFile, "%d-", matrix_size);
    else
        fprintf(pFile, "%s-", obsFPath);
    fprintf(pFile, "%f \n", final_loglik);
    fclose(pFile);

}


void theta_parser2(double *theta_vec, char * kern, int num_params)
    //! parse the theta vector, example: "1:0.5:0.1" -> {1, 0.5, 0.1}
{
    int i = 0;
    if(!strcmp(kern,""))
    {
        for( i = 0 ; i < num_params; i++)
            theta_vec[i] = -1;
    }

    char * token;
    while( (token = strsep(&kern,":")) != NULL )
    {
        if (strcmp(token,"?"))
            theta_vec[i] = strtod(token,NULL);
        else
            theta_vec[i] = -1;
        i++;
    }


}
void write_pred_vector(double *zvec, MLE_data * data, int n, double mserror)

{
    /*!
     * Returns initial_theta, starting_theta, target_theta.
     * @param[in] zvec: measurements vector.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] n: number of spatial locations
     * */

    FILE *pFile;
    pFile = fopen("predict_values.txt","w+");
    if(pFile == NULL) {
        fprintf(stderr,"Cannot access the results path(1)\n");
        exit(0);
    }

    fprintf(pFile, "MSPE: %2.6f\n", mserror);
    int i=0;
    for(i=0;i<n;i++){
        //	printf("zvec[i]: %f\n", zvec[i]);
        fprintf(pFile, "%f\n", zvec[i]);
    }
    fclose(pFile);
}

void write_vectors(double * zvec, MLE_data * data, int n)
    //! store locations, measurements, and log files if log=1
    /*!
     * Returns initial_theta, starting_theta, target_theta.
     * @param[in] zvec: measurements vector.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] n: number of spatial locations 
     * */
{

    int i = 1;
    FILE *pFileZ, *pFileXY;
    location *l = &data->l1;
    struct stat st = {0};
    char * nFileZ  = (char *) malloc(50 * sizeof(char));
    char * temp    = (char *) malloc(50 * sizeof(char));
    char * nFileXY = (char *) malloc(50 * sizeof(char));
    data->nFileLog = (char *) malloc(50 * sizeof(char));	
    //Create New directory if not exist
    if (stat("./synthetic_ds", &st) == -1) 
        mkdir("./synthetic_ds", 0700);

    snprintf(nFileZ, 50, "%s%d%s", "./synthetic_ds/Z_", n,"_");
    snprintf(nFileXY, 50, "%s%d%s", "./synthetic_ds/LOC_", n,"_");
    snprintf(data->nFileLog, 50, "%s%d%s", "./synthetic_ds/log_", n,"_");

    snprintf(temp, 50, "%s%d", data->nFileLog , i);
    while(doesFileExist(temp) == 1)
    {
        i++;
        snprintf(temp, 50, "%s%d", data->nFileLog , i);
    }

    sprintf(temp, "%d", i);
    strcat(nFileZ, temp);
    strcat(nFileXY, temp);
    strcat(data->nFileLog, temp);


    pFileZ = fopen(nFileZ, "w+");
    pFileXY = fopen(nFileXY, "w+");    


    for(i=0;i<n;i++){    
        fprintf(pFileZ, "%f\n", zvec[i]);
        if(l->z ==NULL)
            fprintf(pFileXY, "%f,%f\n", l->x[i], l->y[i]);
        else
            fprintf(pFileXY, "%f,%f,%f\n", l->x[i], l->y[i], l->z[i]);
    }

    fclose(pFileZ);
    fclose(pFileXY);
}

void write_to_thetafile(char * path, double *theta, int num_params,
        int n, double time_per_iter,
        int total_no_iters,
        double prediction_error,
        double mloe,
        double mmom)
    //! write results (only for testing accuracy should be removed in release)
{
    FILE *pFile;
    int i = 0;

    pFile = fopen(path,"a");

    if(pFile == NULL) {
        printf("Cannot access the results path(3)\n");
        exit(0);
    }
    for(i=0; i<num_params; i++)
    {
        fprintf(pFile, "%f,", theta[i]);
    }
    //fprintf(pFile, "%f,", 2*det);
    //fprintf(pFile, "%f,", dotp);
    //fprintf(pFile, "%f,", (double) (n / 2.0) * log(2.0 * PI));
    fprintf(pFile, "%d,", n);
    fprintf(pFile, "%d,", total_no_iters);
    fprintf(pFile, "%f,", time_per_iter);
    fprintf(pFile, "%f,", prediction_error);
    fprintf(pFile, "%f,", mloe);
    fprintf(pFile, "%f\n", mmom);
    fclose(pFile);

}


void write_to_estimatedtheta(char * path, double *theta, int num_params,
        int n, double prediction_time,
        double mloe_mmom_time,
        double prediction_error1,
        double prediction_error2,
        double prediction_error,
        double mloe,
        double mmom,
        int zvecs)
    //! write results (only for testing accuracy should be removed in release)
{
    FILE *pFile;
    int i;
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
        pFile = fopen(path,"a");

        if(pFile == NULL) {
            printf("Cannot access the results path(3)\n");
            exit(0);
        }

        fprintf(pFile, "%d, ", n);
        for(i=0; i<num_params; i++)
        {
            fprintf(pFile, "%f ", theta[i]);
        }
        //fprintf(pFile, "%f,", 2*det);
        //fprintf(pFile, "%f,", dotp);
        //fprintf(pFile, "%f,", (double) (n / 2.0) * log(2.0 * PI));
        fprintf(pFile, "  ---- %f ", prediction_time);
        fprintf(pFile, "%g ", prediction_error1);
        fprintf(pFile, "%g ", prediction_error2);
        fprintf(pFile, "%g --- ", prediction_error);
        fprintf(pFile, "%f ", mloe_mmom_time);       
        fprintf(pFile, "%f ", mloe);
        fprintf(pFile, "%f ", mmom);
        fprintf(pFile, "%d\n", zvecs);
        fclose(pFile);
#if defined(CHAMELEON_USE_MPI)
    }
#endif

}
//void readObsFile(char *obsfile, int n, double * streamdata)
//! Read real observation file
/* @param[in] obsfile: observation file  path.
 * @param[in]  n : number of spatial locations (observations).
 * @param[out] streamdata: observations vector 
 */
/*{

  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  int count = 0;

  fp = fopen(obsfile, "r");
  if (fp == NULL)
  {
  printf("cannot open observations file\n");
  exit(EXIT_FAILURE);
  }


  while ((read = getline(&line, &len, fp)) != -1) 
  streamdata[count++]=atof(line);

  fclose(fp);
  free(line);
  }
  */
void shuffle(double *array, location * locations, size_t n)
    //! shuffle an array
{
    if (n > 1) 
    {
        size_t i;
        //    unsigned int *seed = &exageostat_seed;
        //    size_t jj = i + rand() / (RAND_MAX / (n - i) + 1);
        for (i = 0; i < n - 1; i++) 
        {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);

            double t = array[j];
            array[j] = array[i];
            array[i] = t;
            double xtemp = locations->x[j];
            locations->x[j] = locations->x[i];
            locations->x[i] = xtemp;
            double ytemp = locations->y[j];
            locations->y[j] = locations->y[i];
            locations->y[i] = ytemp;

        }
    }
}

void shuffle2(double *array, double *array2, location * locations, size_t n)
    //! shuffle an array
{
    if (n > 1)
    {
        size_t i;
        //    unsigned int *seed = &exageostat_seed;
        //    size_t jj = i + rand() / (RAND_MAX / (n - i) + 1);
        for (i = 0; i < n - 1; i++)
        {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);

            double t = array[j];
            array[j] = array[i];
            array[i] = t;
            double t2 = array2[j];
            array2[j] = array2[i];
            array2[i] = t2;
            double xtemp = locations->x[j];
            locations->x[j] = locations->x[i];
            locations->x[i] = xtemp;
            double ytemp = locations->y[j];
            locations->y[j] = locations->y[i];
            locations->y[i] = ytemp;

        }
    }
}

void theta_parser(double *initial_theta, double *target_theta, double *starting_theta,
        char *ikernel, char *kernel, double *lb,
        double *up, int test, int num_params)
    //! Parse initial_theta, target_theta, starting_theta, inputs
    /*! 
     * Returns initial_theta, starting_theta, target_theta.
     * @param[out] initial_theta: initial_theta Vector with three parameter (Variance, Range, Smoothness)
     that is used to to generate the Covariance Matrix and initial Z vector.
     * @param[out] starting_theta: theta Vector with three parameter (Variance, Range, Smoothness)
     that is used to to generate the Covariance Matrix of the first MLE iteration.
     * @param[out] target_theta: target theta Vector with three parameter (Variance, Range, Smoothness) unknown theta parameter should be shown as '?'.
     * @param[in] ikernel: initial_theta Vector as string.
     * @param[in] kernel:  target_theta Vector as string.                        
     * @param[in] lb: optimization lower bounds vector ( lb_1, lb_2, lb_3).
     * @param[in] up: optimization upper bounds vector ( ub_1, ub_2, ub_3).
     * @param[in] test: if test=1 ->test running mode, test=0-> real running mode.
     * */
{

    int i = 0;
    //	if (strcmp (initial_theta, "") != 0)
    theta_parser2(initial_theta, ikernel, num_params);
    theta_parser2(target_theta, kernel, num_params);

    for(i = 0; i < num_params; i++)
    {
        if(target_theta[i] != -1)
        {
            lb[i] = target_theta[i];
            up[i] = target_theta[i];
            starting_theta[i] = target_theta[i];
        }
    }

}

void init_optimizer( nlopt_opt * opt, double *lb, double *up, double tol)
    //! Initialize the NLOPT optimizer
    /*!
     * Returns nlopt_opt object.
     * @param[in] lb: optimization lower bounds vector ( lb_1, lb_2, lb_3).
     * @param[in] up: optimization upper bounds vector ( ub_1, ub_2, ub_3).
     * @param[in] tol: a tolerance that is used for the purpose of stopping criteria only.
     * @param[out] opt: nlopt_opt object.
     * */
{
    //initalizing opt library
    nlopt_set_lower_bounds(*opt, lb);
    nlopt_set_upper_bounds(*opt, up);
    nlopt_set_ftol_abs(*opt, tol);
}

void print_summary(int test, int N, int ncores, int gpus, int ts, int lts, char *computation, int zvecs, int p_grid, int q_grid, int precision)
    //! print the summary of MLE inputs.
{
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
        fprintf(stderr,"********************SUMMARY**********************\n");
        if(test == 0)
            fprintf(stderr,"#Execution Mode: Real Dataset\n");
        else
            fprintf(stderr,"#Synthetic Dataset\n");
        fprintf(stderr,"Number of Locations: %d\n", N);
        fprintf(stderr,"#Threads per node: %d\n", ncores);
        fprintf(stderr,"#GPUs: %d\n", gpus);
        if(precision ==0)
            fprintf(stderr,"#Double Precision!\n");
        else if(precision ==1)
            fprintf(stderr,"#Single Precision!\n");
        else if(precision == 2)
            fprintf(stderr,"#Single/Double Precision!\n");
        if (strcmp (computation, "exact") == 0 || strcmp (computation, "diag_approx") == 0)
            fprintf(stderr,"#Dense Tile Size: %d\n", ts);
        else
            fprintf(stderr,"LR Tile Size: %d\n", lts);
        fprintf(stderr,"#%s computation\n", computation);
        fprintf(stderr,"#Obervation Vectors (Z): %d\n", zvecs);
        fprintf(stderr,"p=%d, q=%d\n", p_grid, q_grid);
        fprintf(stderr,"***************************************************\n");
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}


int print_result(MLE_data *data, double *starting_theta, int N, int zvecs, int ncores, int ts, int test, double *initial_theta, char *computation, int p_grid, int q_grid, double final_loglik, double prediction_error)
//! print results (only for testing accuracy should be removed in release)
{
    int num_params;
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif

        if(strcmp(data->kernel_fun, "univariate_matern_stationary")   == 0)
	    num_params = 3; 	
	else if (strcmp(data->kernel_fun, "univariate_pow_exp_stationary")   == 0)
            num_params = 3;
        else if(strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary")   == 0)
            num_params = 4;
        else if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
            num_params = 9;
        else if(strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0)
            num_params = 11;
        else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0)
            num_params = 6;
        else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0)
            num_params = 6;
        else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
            num_params = 6;
        else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
            num_params = 6;
        else if(strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary")   == 0)
            num_params = 7;
        else
        {
            fprintf(stderr,"Choosen kernel is not exist(7)!\n");
            fprintf(stderr, "Called function is: %s\n",__func__);
            exit(0);
        }
        //num_params = strcmp(data->kernel_fun, "stationary_kernel")   == 0? 3 : 9;
        printf("********************SUMMARY*****************************\n");
        printf("#Total Number of Iterations=%d\n",data->iter_count);
        printf("#Total Optimization Time= %6.2f\n", data->total_exec_time);
        printf("#Found Maximum at (");
        int i=0;
        if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        {
            printf("%2.6f, %2.6f,", data->variance1, data->variance2);
            i = 2;
        }
        else
            i = 0;

        for(;i<num_params-1; i++)
            printf("%0.8f, ", starting_theta[i]);
        printf("%0.8f", starting_theta[i]);
        printf(") \n");

        printf("#Number of Locations: %d\n", results.problem_size);
        printf("#Computation: %s\n", results.computation);
        printf("#Kernel: %s\n", results.kernel);
        printf("#Dataset: %s\n", results.ds_type);
        printf("#Precision: %s\n", results.precision);
        printf("#Sample: %d\n", results.z_sample);
        printf("#Dense tile size: %d\n", results.dense_ts);
        printf("#LR tile size (if used): %d\n", results.lr_ts);
        printf("#LR accuracy/max_rank (if used): %d/%d\n", results.lr_acc, results.lr_maxrank);
        printf("#CPU cores per node: %d\n", results.ncores);
        printf("#GPUs per node: %d\n", results.ngpus);
        printf("#p_grid: %d, q_grid: %d\n", results.p, results.q);
        printf("***************************************************\n\n");

        if(data->log == 1)
        {

            fprintf(data->pFileLog,"Total Number of Iterations=%d\n", data->iter_count);
            fprintf(data->pFileLog,"Total Optimization Time= %6.2f secs\n", data->total_exec_time);
            fprintf(data->pFileLog,"Found Maximum at (");
            if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
            {
                printf("%2.6f, %2.6f,", data->variance1, data->variance2);
                i = 2;
            }
            else
                i = 0;

            for(;i<num_params; i++)
                fprintf(data->pFileLog,"%g, ", starting_theta[i]);
            fprintf(stderr,")\n");
        }

        results.time_per_iteration = data->avg_exec_time_per_iter/data->iter_count;
        results.flops_per_iteration     = data->avg_flops_per_iter/data->iter_count;
        results.total_mle_time = data->avg_exec_time_per_iter;

        FILE *pFile;
        pFile = fopen("results.log","a");

        if(pFile == NULL) {
            printf("Cannot access the results path(3)\n");
            return -1;
        }

        fprintf(pFile, "%d ", results.problem_size);
        fprintf(pFile, "%s ", results.computation);
        fprintf(pFile, "%s ", results.kernel);
        //fprintf(pFile, "%s ", results.ds_type);	
        //fprintf(pFile, "%s ", results.precision);
        //fprintf(pFile, "%d ", results.z_sample);
        fprintf(pFile, "%d ", results.dense_ts);
        fprintf(pFile, "%d ", results.lr_ts);
        fprintf(pFile, "%d ", results.lr_acc);
        fprintf(pFile, "%d ", results.lr_maxrank);
        fprintf(pFile, "%d ", results.ncores);
        fprintf(pFile, "%d ", results.ngpus);
        //fprintf(pFile, "%d ", results.p);
        //fprintf(pFile, "%d ", results.q);

        //for(i=0; i<results.num_params; i++)
            //fprintf(pFile, "%6.6f ", results.initial_theta[i]);
        //for(i=0; i<results.num_params; i++)
            //fprintf(pFile, "%6.6f ", results.starting_theta[i]);
        //for(i=0; i<results.num_params; i++)
            //fprintf(pFile, "%6.6f ", results.estimated_theta[i]);
        fprintf(pFile, "%6.6f ", results.final_loglik);
        //fprintf(pFile, "%6.6f ", results.mse_pred1);
        //fprintf(pFile, "%6.6f ", results.mse_pred2);
        fprintf(pFile, "%6.6f ", results.mse_pred);

        fprintf(pFile, "%6.6f ", data->avg_exec_time_gen_stage);
	fprintf(pFile, "%6.6f ", results.time_per_iteration);
	//fprintf(pFile, "%6.6f ", results.total_mle_time);
        fprintf(pFile, "%6.6f ", results.total_pred_time);
	
	fprintf(pFile, "%6.6f ", data->avg_flops_gen_stage);
	fprintf(pFile, "%6.6f ", results.flops_per_iteration);
        fprintf(pFile, "%6.6f ", results.total_pred_flops);

        //fprintf(pFile, "%6.6f ", results.mloe);
        //fprintf(pFile, "%6.2f ", results.mmom);
        //fprintf(pFile, "%s ", results.mloe_exec);
        //fprintf(pFile, "%6.2f ", results.total_mloe_mmom_time);
        //fprintf(pFile, "%6.2f ", results.matrix_gen_mloe_mmom_time);
        //fprintf(pFile, "%6.2f ", results.cho_fact_mloe_mmom_time);
        //fprintf(pFile, "%6.2f ", results.loop_mloe_mmom_time);
        //fprintf(pFile, "%6.2f", results.total_mloe_mmom_flops);
	fprintf(pFile, "\n");



        //		if(zvecs==1)
        //			write_to_file("results.txt", N, ncores, ts, test, computation, data->async, data->obsFPath, data->total_exec_time, data->avg_exec_time_per_iter, data->avg_flops_per_iter, data->precision, data->diag_thick, data->final_loglik, N);



        //		char str[80];
        //		char buf[12];
        //		if(strcmp(data->kernel_fun, "univariate_matern_stationary")   == 0)
        //			strcpy(str, "u-station");
        //		else if(strcmp(data->kernel_fun, "univariate_pow_exp_stationary")   == 0)
        //			strcpy(str, "u-pow-exp-station");
        //		else if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
        //			strcpy(str, "u-non-station");
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0)
        //			strcpy(str, "b-matern-flex");
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0)
        //			strcpy(str, "b-matern-pars");
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0)
        //			strcpy(str, "b-matern-pars2");
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0)
        //			strcpy(str, "b-matern-pars-pro");
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        //			strcpy(str, "b-matern-pars2-pro");

        //		strcat(str, "-");
        //		strcat(str, "theta-");
        //		strcat(str, data->computation);
        //		strcat(str, "-");
        //		if(strcmp (data->computation, "exact") == 0)
        //		{

        //			sprintf(buf, "%d", data->diag_thick);
        //			strcat(str, buf);
        //			strcat(str, "-0-0-");
        //		}
        //		else	if(strcmp (data->computation, "lr_approx") == 0)
        //		{
        //			strcat(str, "0-");
        //			sprintf(buf, "%d", data->hicma_acc);
        //			strcat(str, buf);
        //			strcat(str, "-");
        //			sprintf(buf, "%d", data->hicma_maxrank);
        //			strcat(str, buf);
        //			strcat(str, "-");
        //		}
        //		else if (strcmp (data->computation, "diag_approx") == 0)
        //		{
        //			sprintf(buf, "%d", data->diag_thick);
        //			strcat(str, buf);
        //			strcat(str, "-0-0-");
        //		}

        //		if(data->async==0)
        //			strcat(str, "syn-");
        //		else
        //			strcat(str, "asyn-");
        //		if(data->precision == 0)
        //			strcat(str, "d_prec-");
        //		else if(data->precision ==1)
        //			strcat(str, "s_prec-");
        //		else if(data->precision == 2)
        //		{
        //			strcat(str, "sd_prec-");
        //		}
        //
        //		if(strcmp(data->locsFPath,"") == 0)
        //		{
        //			for(int i=0; i<6; i++)
        //			{

        //				sprintf(buf, "%0.3f-", initial_theta[i]);
        //				strcat(str, buf);
        //			}
        //		}
        //		else
        //		{		
        //			strcat(str, "real-");
        //		}

        //		sprintf(buf, "%d", N);
        //		strcat(str, buf);
        //		strcat(str, ".txt");

        //		printf("%s\n", str);
        //exit(0);
        //		if(strcmp(data->kernel_fun, "univariate_matern_stationary")   == 0 || strcmp(data->kernel_fun, "univariate_pow_exp_stationary")   == 0)
        //			num_params = 3;
        //		else if(strcmp(data->kernel_fun, "univariate_matern_non_stationary")   == 0)
        //			num_params = 9;
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_flexible")   == 0)
        //			num_params = 11;
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious")   == 0)
        //			num_params = 6;
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious2")   == 0)
        //			num_params = 6;
        //		else if(strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile")   == 0 || strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile")   == 0)
        //		{
        //			starting_theta[0] = data->variance1;
        //			starting_theta[1] = data->variance2;
        //			num_params = 6;
        //		}
        //		else
        //		{ 
        //			fprintf(stderr,"Choosen kernel is not exist(5)!\n");
        //			fprintf(stderr, "Called function is: %s\n",__func__);
        //			exit(0);
        //		}
        //num_params = strcmp(data->kernel_fun, "stationary_kernel")   == 0? 3 : 9;
        //		write_to_thetafile(str, starting_theta, num_params, zvecs, data->avg_exec_time_per_iter,data->iter_count, prediction_error, data->mloe, data->mmom );
        fclose(pFile);
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}


double cWtime(void)
    //! get time
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

//void readlocfile(char* loc_file, int n,  location* l1)
//! Read real location file
/*!
 * Returns locations l1.
 * @param[in] loc_file: locations file  path.
 * @param[in]  n : number of spatial locations (observations).
 * @param[out] l1 :  location struct with (x,y) locations.
 */
/*{
  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  char * pch;
  int i=0;
  fp = fopen(loc_file, "r");
  if (fp == NULL)
  {
  printf("cannot read locations file: %s\n", loc_file);
  exit(EXIT_FAILURE);
  }
  while ((read = getline(&line, &len, fp)) != -1 && i < n) {
  pch = strtok (line,",");
  while (pch != NULL)
  {
  l1->x[i] = atof(pch);
  pch = strtok (NULL, ",");
  l1->y[i] = atof(pch);
  pch = strtok (NULL, ",");
  }
  i++;
  }

  fclose(fp);
  }
  */
void write_prediction_result(char * path, int matrix_size,int no_missing, double MSE1, double MSE2, double MSE, double solve_time, double flops)
    //! write prediction results (only for testing accuracy should be removed in release).
{
    FILE *pFile;

    pFile = fopen(path,"a");

    if(pFile == NULL) {

        fprintf(stderr,"Cannot access the results path(2)\n");
        exit(0);
    }

    fprintf(pFile, "%d\t", matrix_size);
    fprintf(pFile, "%d\t", no_missing);
    fprintf(pFile, "%f\t", MSE1);
    fprintf(pFile, "%f\t", MSE2);
    fprintf(pFile, "%f\t", MSE);
    fprintf(pFile, "%6.2f secs\t", solve_time);
    fprintf(pFile, "%6.2f Gflops\n", flops);	

    fclose(pFile);

}

// check if the file exist or not
int doesFileExist(const char *filename)
    /*! variables from argument
     * Returns 0 if file not exist  and 1 if file exist.
     * @param[in] filename: file path.
     */
{
    struct stat st;
    int result = stat(filename, &st);
    return result == 0;
}

void init_log (MLE_data * data)
    //! init and open log files if log==1
    /*! 
     * Returns MLE_data struct with new log files.
     * @param[out] data: MLE_data struct with different MLE inputs.
     */
{
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
        data->pFileLog = fopen(data->nFileLog, "w+");
        fprintf (data->pFileLog, "\t\tlog file is generated by ExaGeoStat application\n"); 
        fprintf (data->pFileLog, "\t\t============================================\n");
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}

void finalize_log (MLE_data * data)
    //! finalize and close log files if log==1
    /*!
     * Returns MLE_data struct with new log files.
     * @param[out] data: MLE_data struct with different MLE inputs.
     */
{
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
        fclose(data->pFileLog);
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}

/*
//check accuracy
acc_struct check_acc(MLE_data *data, int n, int ts)
{
#if defined(EXAGEOSTAT_USE_HICMA)
int main_print_index    = 0;
int print_mat           = 0;
int M                   = n;
int N                   = n /ts * data->hicma_maxrank;
int MB                  = ts;
int main_print_mat      = 0;
int LDA                 = M;
int check               = 1;
int set_diag            = 0;
double diagVal          = M;

acc_struct acc;
MORSE_desc_t *descAD    = data->hicma_descCD;
MORSE_desc_t *descAUV   = data->hicma_descCUV;
MORSE_desc_t *descDense = data->descC;
MORSE_desc_t *descArk   = data->hicma_descCrk;
double *Adense2 = data->Adense2;
double *Adense  = data->Adense;
double *Ahicma = (double *) malloc(n * n * sizeof(double));
double *AhicmaT = (double *) malloc(n * n * sizeof(double));


HICMA_zuncompress(MorseLower, descAUV, descDense, descArk);
HICMA_zdiag_vec2mat(descAD, descDense);
//PASTE_CODE_FREE_MATRIX( descAD ); //@KADIRLBL001
descAD = descDense; // descAD was only diagonals.
// After this line, descAD is dense matrix containing approximate L
// So no need to adapt below code for descAD containg only diagonals.

////PROGRESS("checking accuracy");
if( MORSE_My_Mpi_Rank()==0){
#ifndef COMPLEX
double one = 1.0, zero = 0.0, minusone = -1.0, diagVal = M;
if(main_print_mat){printf("data->Adense2\n");/*printmat(data->Adense2,M,M,LDA,MB, MB);*/ /*}
                                                                                           double* swork = calloc(2*M, sizeof(double));
                                                                                           double normA;
                                                                                           {size_t i, j;
                                                                                           for(j = 0; j < M; j++){
                                                                                           for(i = 0; i < j; i++){
                                                                                           Adense2[j*LDA+i] = zero;
                                                                                           }
                                                                                           }
                                                                                           }

//PROGRESS("normaA started");
HICMA_znormest(M, M, data->Adense2, &normA, swork);
// Ahicma: result of TLR potrf
MORSE_Tile_to_Lapack(descAD, Ahicma, n);
//PASTE_TILE_TO_LAPACK( descAD, Ahicma, check, double, LDA, M );
{size_t i, j;
for(j = 0; j < M; j++){
for(i = 0; i < j; i++){
Ahicma[j*LDA+i] = zero;
}
}
}

if(set_diag){size_t j; for(j = 0; j < M; j++){ Ahicma[j*LDA+j] = diagVal; } }
//if(main_print_mat){printf("Ahicma\n");printmat(Ahicma,M,M,LDA, MB, MB);}

//LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', M, Ahicma, LDA);
// AhicmaT: transpose of Ahicma
//PROGRESS("copy descAd into AhicmaT started");

MORSE_Tile_to_Lapack(descAD, AhicmaT, n);
//PASTE_TILE_TO_LAPACK( descAD,  AhicmaT, check, double, LDA, M );

{size_t i, j;
for(j = 0; j < M; j++){
for(i = 0; i < j; i++){
data->Adense[j*LDA+i] = zero;
}
}
}


//if(main_print_mat){printf("Ahicma-upperzero\n");printmat(Ahicma,M,M,LDA, MB, MB);}
//PROGRESS("Transpose A started");
//LAPACKE_dge_trans(LAPACK_COL_MAJOR, M, M, Ahicma, LDA, AhicmaT, LDA);
//if(main_print_mat){printf("AhicmaT\n");printmat(AhicmaT,M,M,LDA, MB, MB);}
//PROGRESS("TRMM started");
cblas_dtrmm (CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, M, one, Ahicma, LDA, AhicmaT, LDA);
//if(main_print_mat){printf("Ahicma*AhicmaT\n");printmat(AhicmaT,M,M,LDA, MB, MB);}
//double tmpnorm;normest(M, M, AhicmaT, &tmpnorm, swork);printf("tmpnorm:%e\n",tmpnorm);
{size_t i, j;
for(j = 0; j < M; j++){
for(i = 0; i < j; i++){
AhicmaT[j*LDA+i] = zero;
}
}
}

size_t nelm = M * M;
if(main_print_mat)printf("nelm:%zu M:%d N:%d\n", nelm, M, N);
//PROGRESS("DAXPY started");
cblas_daxpy(nelm, minusone, AhicmaT, 1, data->Adense, 1);
//if(main_print_mat){printf("data->Adense-(Ahicma*AhicmaT)\n");printmat(data->Adense,M,M,LDA, MB, MB);}


//PROGRESS("Norm of difference started");
HICMA_znormest(M, M, data->Adense, &acc.normDenseAppDiff, swork);
acc.accuracyDenseAppDiff = acc.normDenseAppDiff/normA;
acc.normA = normA;
//acc.normDenseAppDiff = normDenseAppDiff;
//printf("normA:%.2e normDenseAppdiff:%.2e Accuracy: %.2e\n", normA, normDenseAppDiff,  accuracyDenseAppDiff);
#endif
} else {

    MORSE_Tile_to_Lapack(descAD, Ahicma, n);
    MORSE_Tile_to_Lapack(descAD, AhicmaT, n);

    //PASTE_TILE_TO_LAPACK( descAD, Ahicma, check, double, LDA, M );
    //PASTE_TILE_TO_LAPACK( descAD,  AhicmaT, check, double, LDA, M );
}
//PROGRESS("checking accuracy is finished");


//  PASTE_CODE_FREE_MATRIX( descAUV );
//PROGRESS("descAUV is freed");
//if(check == 0) { // If there is no check, then descAD and descDense are different. Refer to @KADIRLBL001
//  PASTE_CODE_FREE_MATRIX( descAD );
//PROGRESS("descAD is freed");
// }
//    PASTE_CODE_FREE_MATRIX( descArk );
//PROGRESS("descArk is freed");

// PASTE_CODE_FREE_MATRIX( descDense );
//PROGRESS("descDense is freed");
//PROGRESS("freed descs");
return acc;
#endif
}
*/

void split_data(MLE_data * data, location *locations, double *Z, double *Zactual,  int *N, int nZmiss)
{

    int i = 0;
    int nZobs = *N-nZmiss;
    shuffle(Z, locations, *N);
    data->l1.x     = (double *) malloc(nZobs * sizeof(double));
    data->l1.y     = (double *) malloc(nZobs * sizeof(double));
    data->lmiss.x     = (double *) malloc(nZmiss * sizeof(double));
    data->lmiss.y     = (double *) malloc(nZmiss * sizeof(double));

    for( i = 0; i < nZobs; i++)
    {
        data->l1.x[i] = locations->x[i];
        data->l1.y[i] = locations->y[i];
    }

    data->lobs.x = data->l1.x;
    data->lobs.y = data->l1.y;

    for( i = 0; i < nZmiss; i++)
    {
        data->lmiss.x[i] = locations->x[i + nZobs];
        data->lmiss.y[i] = locations->y[i + nZobs];
        Zactual[i] = Z[i + nZobs];
    }
    *N = nZobs;

}

void pick_random_points_noshuffle(MLE_data *data, double *Zobs, double *Zactual, int nZmiss, int nZobs, int N)
{
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double *Z;
    int i = 0;

    //memory allocation
    Z       = (double *) malloc(N * sizeof(double));
    l.x     = (double *) malloc(N/2 * sizeof(double));
    l.y     = (double *) malloc(N/2 * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    for( i = 0; i < N/2 ; i++)
    {
        l.x[i]=data->l1.x[i];
        l.y[i]=data->l1.y[i];
    }

    //        shuffle(Z, &l, N);
    for( i = 0; i < nZobs ; i++)
        Zobs[i] = Z[nZmiss+i];

    for ( i = 0; i < nZmiss; i++)
        Zactual[i] = Z[i];

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x  = &l.x[nZmiss];
    lobs->y  = &l.y[nZmiss];

    locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    free(Z);
}

void pick_random_points(MLE_data *data, double *Zobs, double *Zactual, int nZmiss, int nZobs, int N)
{
    //initialization
    location l;
    location *lmiss;
    location *lobs;


    double *Z;
    int i = 0;

    //memory allocation
    Z       = (double *) malloc(N * sizeof(double));
    l.x     = (double *) malloc(N * sizeof(double));
    l.y     = (double *) malloc(N * sizeof(double));


    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    for( i = 0; i < N ; i++)
    {
        l.x[i]=data->l1.x[i];
        l.y[i]=data->l1.y[i];
    }

    shuffle(Z, &l, N);
    for( i = 0; i < nZobs ; i++)
        Zobs[i] = Z[nZmiss+i];

    for ( i = 0; i < nZmiss; i++)
        Zactual[i] = Z[i];



    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x  = &l.x[nZmiss];
    lobs->y  = &l.y[nZmiss];

    locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    free(Z);
}


void pick_random_points2(MLE_data *data, double *Zobs, double *Zactual, int nZmiss, int nZobs, int N)
{
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double *Z;
    double *Z1;
    double *Z2;
    int i = 0;
    int p = 2;




    //memory allocation
    Z       = (double *) malloc(N * sizeof(double));
    l.x     = (double *) malloc(N/p * sizeof(double));
    l.y     = (double *) malloc(N/p * sizeof(double));
    Z1       = (double *) malloc(N * sizeof(double));
    Z2       = (double *) malloc(N * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    double sum=0;
    int j=0;
    for(i=0;i<N;i+=2)
    {
        Z1[j]=Z[i];
        Z2[j++]=Z[i+1];
    }

    /*
       for(i=0;i<N;i++)
       {
       printf("%d - %f - %f - %f\n", i, Z[i], Z1[i], Z2[i]);
       }

       exit(0);
       */
    for( i = 0; i < N/p ; i++)
    {
        l.x[i]=data->l1.x[i];
        l.y[i]=data->l1.y[i];
    }

    shuffle2(Z1, Z2, &l, N/p);



    //actual vector
    j=0;
    for( i = 0; i < nZmiss ; i++)
    {
        Zactual[j] =Z1[i];
        Zactual[j+1]=Z2[i];
        j+=2;
    }

    //	printf("nZmiss: %d, N:%d\n", nZmiss, N);
    //observation vector
    j=0;
    for( i = 0; i < nZobs ; i++)
    {
        Zobs[j]=Z1[nZmiss+i];
        Zobs[j+1]=Z2[nZmiss+i];
        j+=2;
    }

    /*	printf("nZobs=%d\n", nZobs);
        for(i=0;i<2*nZmiss+1;i++)
        {
        printf("%f - %f - %f\n", Z1[i], Z2[i], Zactual[i]);
        }

        exit(0);
        */

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x  = &l.x[nZmiss];
    lobs->y  = &l.y[nZmiss];

    //TODO: check why using them give wrong answer
    //	locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    //	locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    //	locations_obs_zsort_inplace(nZobs, lobs, &Zobs[nZobs/2]);
    //	locations_obs_zsort_inplace(nZmiss, lmiss, &Zactual[nZmiss/2]);
    //	free(Z);
    free(Z1);        free(Z2);
}

void generate_interior_points(MLE_data *data, double *Zobs, double *Zactual, int nZmiss, int nZobs, int N)
{
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double *Z;
    int i = 0;

    if(nZmiss >= nZobs)
    {
        fprintf(stderr,"Cannot generate missing locations larger than or equal the observed ones\n");
        return;
    }

    //memory allocation
    Z       = (double *) malloc(N * sizeof(double));
    l.x     = (double *) malloc(N * sizeof(double));
    l.y     = (double *) malloc(N * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);


    for( i = 0; i < N ; i++)
    {
        l.x[i]=data->l1.x[i];
        l.y[i]=data->l1.y[i];
    }

    shuffle(Z, &l, N);

    for( i = 0; i < nZobs ; i++)
        Zobs[i] = Z[nZmiss+i];


    lmiss = &(data->lmiss);
    lobs  = &(data->lobs);

    for (i =0; i<nZmiss; i++)
    {
        l.x[i] = (l.x[i] + l.x[i+1]) / 2.0;
        l.y[i] = (l.y[i] + l.y[i+1]) / 2.0;
    }

    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x  = &l.x[nZmiss];
    lobs->y  = &l.y[nZmiss];
    free(Z);
}

// This function converts decimal degrees to radians
static double deg2rad(double deg) {
    return (deg * PI / 180);
}
//  This function converts radians to decimal degrees
static double rad2deg(double rad) {
    return (rad * 180 / PI);
}
/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
static double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r)/2);
    v = sin((lon2r - lon1r)/2);
    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

static double calculateDistance(double x1, double y1, double x2, double y2, int distance_metric) {

    if(distance_metric == 1)
        return distanceEarth(x1, y1, x2, y2);
    return  sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}




double core_matern_vector (double x0, double y0, double x1, double y1, double *localtheta, int distance_metric) {

    int i, j;
    double expr = 0.0;
    double con = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2,(localtheta[2]-1)) * tgamma(localtheta[2]);
    con = 1.0/con;
    con = sigma_square * con;

    expr = calculateDistance(x0, y0, x1, y1, distance_metric)/localtheta[1];
    //printf("expr = %f\n", expr);
    if(expr == 0)
        return sigma_square /*+ 1e-4*/;
    else
        return con*pow(expr, localtheta[2])*gsl_sf_bessel_Knu(localtheta[2],expr); // Matern Function
}

void fwrite_array(int m, int n, int ld, double* arr, char* file){
    FILE* fp = fopen(file, "w");
    if(fp == NULL){
        fprintf(stderr, "File %s cannot be opened to write\n", file);
        exit(1);
    }
    int i, j;
    fprintf(fp, "%d %d\n", m, n);
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            fprintf(fp, "%d\t", (int)arr[ld*j+i] );
        }
        fprintf(fp, "\n" );
    }
    fclose(fp);
}

void print_array(int m, int n, int ld, double* arr, FILE* fp){
    int i, j;
    fprintf(fp, "%d %d\n", m, n);
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            fprintf(fp, "%d\t", (int)arr[ld*j+i] );
        }
        fprintf(fp, "\n" );
    }
}
void checkpointing(char *path, int iter_count, double* theta, double loglik, int num_params)
    //! write prediction results (only for testing accuracy should be removed in release).
{
    FILE *pFile;
    int i;

    pFile = fopen(path,"a");

    if(pFile == NULL) {

        fprintf(stderr,"Cannot access the results path\n");
        exit(0);
    }

    fprintf(pFile, "%d ", iter_count);
    for(i=0; i<num_params; i++)
    {
        fprintf(pFile, "%.17g ", theta[i]);
    }
    fprintf(pFile, "%.17g ", loglik);

    fprintf(pFile, "\n");
    fclose(pFile);

}
bool recover(char *path, int iter_count, double* theta, double* loglik, int num_params)
{

    FILE * fp;
    char * line     = NULL;
    size_t len      = 0;
    ssize_t read;
    int count       = 0;
    int i = 0;
    char *pch;
    fp = fopen(path, "r");
    if (fp == NULL)
    {
        printf("cannot open observations file\n");
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &len, fp)) != -1) {
        pch = strtok(line, " ");
        count =atoi(pch);
        if(count == iter_count)
        {
            pch = strtok (NULL, " ");
            for(i=0;i<num_params;i++)
            {
                theta[i] = atof(pch);
                pch = strtok (NULL, " ");
            }
            *loglik = atof(pch);
            fclose(fp);
            free(line);
            //			printf("recover: %d %f %f %f\n", count, theta[0], theta[1], theta[2], *loglik);
            return true;

        }
        count++;
    }



    fclose(fp);
    free(line);

    return false;

}


double* pred_idw(MLE_data *data, double *z_miss, double *z_actual, double*z_obs, int nZmiss, int nZobs)
{
    int i = 0, j = 0;
    double sigma_1 = 0;
    double sigma_2 = 0;
    double error = 0, dij = 0;
    double error1 = 0, error2 = 0;
    double x1, y1, x2, y2;
    int n = 2;
    location* lmiss= &data->lmiss;
    location* lobs= &data->lobs;
    double* mspe     = (double *) malloc(3 * sizeof(double));
    for (j=0; j<nZmiss; j++)
    {
        x2=lmiss->x[j];
        y2=lmiss->y[j];
        for(i=0; i<nZobs; i++)
        {
            x1=lobs->x[i];
            y1=lobs->y[i];
            dij = calculateDistance( x1, y1, x2, y2, 0);
            if(dij != 0)
            {
                sigma_1 += z_obs[i]/pow(dij, n);
                sigma_2 += 1.0 / pow(dij, n);
            }
        }
        z_miss[j] = sigma_1/sigma_2;
        if(j % 2 == 0 )
            error1 += pow((z_actual[j]-z_miss[j]), 2);
        else
            error2 += pow((z_actual[j]-z_miss[j]), 2);
        error += pow((z_actual[j]-z_miss[j]), 2);
        sigma_1 = 0;
        sigma_2 = 0;
    }

    mspe[0]=error/nZmiss;
    mspe[1]=error1/(nZmiss/2);
    mspe[2]=error2/(nZmiss/2);
    return mspe;
}

