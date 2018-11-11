/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
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
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
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

	data->computation 		= "exact";	///< exact or approx computation.
        data->test   		        = 0;      	///< exact or approx computation.
        data->iter_count        	= 0;  		///< number of iterations to converge.
        data->l1.x                      = NULL;   	///< 2D locations for the first dataset (X vector).
	data->l1.y			= NULL;		//< 2D locations for the first dataset (Y vector).
        data->lmiss.x           	= NULL;		///< 2D locations for the missing data (X vector) (prediction stage).
        data->lmiss.y             	= NULL;		///< 2D locations for the missing data (Y vector) (prediction stage).
        data->lobs.x            	= NULL;		///< 2D locations for the observed data (X vector) (prediction stage).
        data->lobs.y             	= NULL;		///< 2D locations for the observed data (Y vector) (prediction stage).
	data->descC               	= NULL;		///< Covariance matrix C descriptor.
        data->descZ             	= NULL;		///< Measurements Z descriptor.
        data->Adense             	= NULL;		///< Dense matrix descriptor in the case of approximation mode - accuracy check.
        data->Adense2              	= NULL;		///< Dense matrix descriptor2 in the case of approximation mode - accuracy check.
        data->descZcpy                  = NULL;		///< A copy of Measurements Z descriptor.
        data->descdet                 	= NULL;		///< Determinant descriptor.
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
        data->actualZFPath              = "";		///< Actual observations file path -- in the case of prediction.
        data->actualZLocFPath           = "";		///< Actial locations file path -- in the case of prediction.
        data->det                       = 0.0;		///< determinant value.
        data->dotp                      = 0.0;		///< dot product value.
        data->mserror                   = 0.0;		///< Mean Square Error (MSE) value.
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
	data->opt_tol			= 5;            ///< The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
        data->opt_max_iters		= -1;     	///< Maximum number of mle iterations.
        data->ooc                       = 0;            ///< Support Out-Of-Core execution, 0-->no, 1-->yes.
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
        //printf("%f %f %f -> %u %u\n", (double)UINT16_MAX, locations->x[i], locations->y[i], x, y);
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
//	printf("%f %f ", locations->x[i], locations->y[i]);
        locations->x[i] = (double)x/(double)UINT16_MAX;
        locations->y[i] = (double)y/(double)UINT16_MAX;
        z[i]		= z_struct[i].z;
//	printf("%f %f %f\n", locations->x[i], locations->y[i], y, z[i]);
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
        int i = 0 ,index = 0, j = 0;
       // unsigned int *seed = &exageostat_seed;
        srand(seed);
	location* locations = (location *) malloc( sizeof(location*));
        //Allocate memory
        locations->x            = (double *) malloc(n * sizeof(double));
        locations->y            = (double *) malloc(n * sizeof(double));
       // if(strcmp(locs_file, "") == 0)
       // {

        	int sqrtn = sqrt(n);        

        	//Check if the input is square number or not
        	if(pow(sqrtn,2) != n)    
                {
                	printf("Please use a perfect square number to generate a valid synthetic dataset.....\n\n");
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
                                locations->y[index++] = (grid[j]-0.5+uniform_distribution(-0.4, 0.4))/sqrtn;
                        }
                free(grid);
		zsort_locations(n, locations);
		return locations;
      	// } 
	/*else
    	{

	        FILE * fp;
       		char * line = NULL;
        	size_t len = 0;
        	ssize_t read;
        	char * pch;

        	fp = fopen(locs_file, "r");
        	if (fp == NULL)
        	{
            		printf("cannot read locations file\n");
            		exit(EXIT_FAILURE);
        	}
        
		while ((read = getline(&line, &len, fp)) != -1) {
            		pch = strtok(line, ",");
            		while (pch != NULL)
            		{
                		locations->x[i] = atof(pch);
                		pch = strtok (NULL, ",");
                		locations->y[i] = atof(pch);
               		 	pch = strtok (NULL, ",");
            		}
			i++;
       		}
	

        	fclose(fp);
        	if (line)
            		free(line);

    	}

	return 0;
*/
}

void print_matrix(char* desc, int m, int n, double* a, int lda)
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


void write_to_file(char * path, int matrix_size,int ncores,int tile_size, int test, char * ikernel, char *computation, int async, char *obsFPath,double total_exec_time,double avg_exec_time_per_iter, double avg_flops_per_iter , int p_grid, int q_grid, double final_loglik, int n)
//! write results in  detail (only for testing accuracy should be removed in release)
{
	FILE *pFile;
	double peakperfomrance;
	double percent;    
	pFile = fopen(path,"a");
	if(pFile == NULL) {
		fprintf(stderr,"Cannot access the results path\n");
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
	fprintf(pFile, "%s-", ikernel);
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


void theta_parser2(double * theta_vec,char * kern)
//! parse the theta vector, example: "1:0.5:0.1" -> {1, 0.5, 0.1}
{
	int i = 0;
	if(!strcmp(kern,""))
    	{
		for( i = 0 ; i < 3 ; i++)
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
        snprintf(nFileXY, 50, "%s%d%s", "./synthetic_ds/XY_", n,"_");
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
                fprintf(pFileXY, "%f,%f\n", l->x[i], l->y[i]);
	}
       
     	fclose(pFileZ);
        fclose(pFileXY);
}

void write_to_thetafile(char * path, double theta0,double theta1, double theta2, double loglik, int n)
//! write results (only for testing accuracy should be removed in release)
{
        FILE *pFile;

        pFile = fopen(path,"a");

        if(pFile == NULL) {
                printf("Cannot access the results path\n");
                exit(0);
        }
        fprintf(pFile, "%f,", theta0);
        fprintf(pFile, "%f,", theta1);
        fprintf(pFile, "%f,", theta2);
        //fprintf(pFile, "%f,", 2*det);
        //fprintf(pFile, "%f,", dotp);
        //fprintf(pFile, "%f,", (double) (n / 2.0) * log(2.0 * PI));
	fprintf(pFile, "%d,", n);
	fprintf(pFile, "%f\n", loglik);
        fclose(pFile);

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
//	unsigned int *seed = &exageostat_seed;
//	size_t jj = i + rand() / (RAND_MAX / (n - i) + 1);
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



void theta_parser(double *initial_theta, double *target_theta, double *starting_theta, char *ikernel, char *kernel, double *lb, double *up, int test)
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
    if (test == 1)
                theta_parser2(initial_theta, ikernel);
        theta_parser2(target_theta, kernel);

        for(i = 0; i < 3; i++)
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
        *opt=nlopt_create( NLOPT_LN_BOBYQA, 3);
        nlopt_set_lower_bounds(*opt, lb);
        nlopt_set_upper_bounds(*opt, up);
        nlopt_set_xtol_rel(*opt, tol);
}

void print_summary(int test, int N, int ncores, int gpus, int ts, char *computation, int zvecs, int p_grid, int q_grid)
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
		fprintf(stderr,"Tile Size: %d\n", ts);
		fprintf(stderr,"#%s computation\n", computation);
		fprintf(stderr,"#Obervation Vectors (Z): %d\n", zvecs);
		fprintf(stderr,"p=%d,q=%d\n", p_grid, q_grid);
		fprintf(stderr,"***************************************************\n");
	#if defined(CHAMELEON_USE_MPI)
	}
	#endif
}


void print_result(MLE_data *data, double *starting_theta, int N, int zvecs, int ncores, int ts, int test, char * ikernel, char *computation, int p_grid, int q_grid, double final_loglik)
//! print results (only for testing accuracy should be removed in release)
{
	#if defined(CHAMELEON_USE_MPI)
	if ( MORSE_My_Mpi_Rank() == 0 )
   	 {
    	#endif
		fprintf(stderr,"Total Number of Iterations=%d\n",data->iter_count);
                fprintf(stderr,"Total Optimization Time= %6.2f\n", data->total_exec_time);
                fprintf(stderr,"Found Maximum at f(%g, %g, %g) \n", starting_theta[0], starting_theta[1], starting_theta[2]);

		if(data->log == 1)
		{
			fprintf(data->pFileLog,"Total Number of Iterations=%d\n", data->iter_count);
			fprintf(data->pFileLog,"Total Optimization Time= %6.2f\n", data->total_exec_time);
			fprintf(data->pFileLog,"Found Maximum at f(%g, %g, %g) \n", starting_theta[0], starting_theta[1], starting_theta[2]);	
		}

                data->avg_exec_time_per_iter = data->avg_exec_time_per_iter/data->iter_count;
                data->avg_flops_per_iter     = data->avg_flops_per_iter/data->iter_count;

            if(zvecs==1)
                    write_to_file("results.txt", N, ncores, ts, test, ikernel, computation, data->async, data->obsFPath, data->total_exec_time, data->avg_exec_time_per_iter, data->avg_flops_per_iter, p_grid, q_grid, data->final_loglik, N);

                write_to_thetafile("theta.txt", starting_theta[0], starting_theta[1], starting_theta[2], data->hicma_acc, N);
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
void write_prediction_result(char * path, int matrix_size,int no_missing, double MSE, double solve_time, double flops)
//! write prediction results (only for testing accuracy should be removed in release).
{
        FILE *pFile;

        pFile = fopen(path,"a");

        if(pFile == NULL) {

                fprintf(stderr,"Cannot access the results path\n");
                exit(0);
        }

        fprintf(pFile, "%d\t", matrix_size);
        fprintf(pFile, "%d\t", no_missing);
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
//                        locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

			free(Z);
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
