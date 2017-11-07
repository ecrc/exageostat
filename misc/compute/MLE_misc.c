/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
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
 * @version 0.1.0
 *
 * @author Sameh Abdulah
 * @date 2017-11-07
 *
 **/
#include "../include/MLE_misc.h"
//***************************************************************************************

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

double uniform_distribution(double rangeLow, double rangeHigh) 
//! Generate uniform distribution between rangeLow , rangeHigh
{
    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = rangeHigh - rangeLow;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

int GenerateXYLoc(int n, char * locs_file, location * locations) 
//! Generate XY location for exact computation (MOORSE)        
{
	//initalization
        int i = 0 ,index = 0, j = 0;
        srand(0);
        if(strcmp(locs_file, "") == 0)
        {

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
		zsort_locations(n,locations);
      	} 
	else
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
	fprintf(pFile, "%d,", n);
        fprintf(pFile, "%f,", theta0);
        fprintf(pFile, "%f,", theta1);
        fprintf(pFile, "%f,", theta2);
        //fprintf(pFile, "%f,", 2*det);
        //fprintf(pFile, "%f,", dotp);
        //fprintf(pFile, "%f,", (double) (n / 2.0) * log(2.0 * PI));
	fprintf(pFile, "%f\n", loglik);
        fclose(pFile);

}

void readObsFile(char *obsfile, int n, double * streamdata)
//! Read real observation file
/* @param[in] obsfile: observation file  path.
* @param[in]  n : number of spatial locations (observations).
* @param[out] streamdata: observations vector 
*/
{

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

void shuffle(double *array, location * locations, size_t n)
//! shuffle an array
{
    if (n > 1) 
    {
        size_t i;
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

                write_to_thetafile("theta.txt", starting_theta[0], starting_theta[1], starting_theta[2], data->final_loglik, N);
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

void readlocfile(char* loc_file, int n,  location* l1)
//! Read real location file
/*!
* Returns locations l1.
* @param[in] loc_file: locations file  path.
* @param[in]  n : number of spatial locations (observations).
* @param[out] l1 :  location struct with (x,y) locations.
*/
{
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



