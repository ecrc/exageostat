/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file flat_file.c
 *
 * Auxiliary functions that are used to read and process flat files.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/flat_file.h"
//***************************************************************************************
location* readLocsFile(char *locs_file, int n)
    //! Read 2D locations from a given
    /*! flatfile
     * Returns location struct.
     * @param[in] locs_file: 2D location file.
     * @param[in]  n : number of spatial locations.
     * */
{

    FILE *fp;
    int i = 0;
    char *line = NULL;
    size_t len  = 0;
    ssize_t read;
    char *pch;        
    location *locations; 


    fp = fopen(locs_file, "r");
    if (fp == NULL)
    {
        printf("cannot read locations file\n");
        printf("%s: \n",locs_file);
        return NULL;
    }
    else
    {
        //Allocate memory
        locations		= (location *) malloc(sizeof(location*));
        locations->x            = (double *) malloc(n * sizeof(double));
        locations->y            = (double *) malloc(n * sizeof(double));
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
    //    zsort_locations(n,locations);
    return locations;
}



location* readLocsFile3d(char *locs_file, int n)
    //! Read 2D locations from a given
    /*! flatfile
     * Returns location struct.
     * @param[in] locs_file: 2D location file.
     * @param[in]  n : number of spatial locations.
     * */
{

    FILE *fp;
    int i = 0;
    char *line = NULL;
    size_t len  = 0;
    ssize_t read;
    char *pch;
    location *locations;


    fp = fopen(locs_file, "r");
    if (fp == NULL)
    {
        printf("cannot open location file%s\n", locs_file);
        exit(EXIT_FAILURE);
    }
    else
    {
        //Allocate memory
        locations               = (location *) malloc(sizeof(location*));
        locations->x            = (double *) malloc(n * sizeof(double));
        locations->y            = (double *) malloc(n * sizeof(double));
        locations->z            = (double *) malloc(n * sizeof(double));      
    }

    while ((read = getline(&line, &len, fp)) != -1) {
        pch = strtok(line, ",");
        while (pch != NULL)
        {
            locations->x[i] = atof(pch);
            pch = strtok (NULL, ",");
            locations->y[i] = atof(pch);
            pch = strtok (NULL, ",");
            locations->z[i] = atof(pch);
            pch = strtok (NULL, ",");
        }
        i++;
    }
    fclose(fp);
    if (line)
        free(line);
    //    zsort_locations(n,locations);
    return locations;
}
double * readObsFile(char *obsfile, int n)
    //! Read real observation file
    /* Returns observations struct.
     *  @param[in] obsfile: observation file  path.
     *  @param[in]  n : number of spatial locations (observations).
     *  @param[out] streamdata: observations vector.
     */
{

    FILE * fp;
    char * line     = NULL;
    size_t len      = 0;
    ssize_t read;
    int count       = 0;

    double *z_vec = (double *) malloc(n * sizeof(double));

    fp = fopen(obsfile, "r");
    if (fp == NULL)
    {
        printf("readObsFile:cannot open observations file: %s\n", obsfile);
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &len, fp)) != -1)
        z_vec[count++]=atof(line);

    fclose(fp);
    free(line);

    return z_vec;
}

