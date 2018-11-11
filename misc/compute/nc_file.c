/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file nc_file.c
 *
 * Auxiliary functions that are used to read and process netCDF files.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#include "../include/nc_file.h"
//***************************************************************************************
int openFileNC(MLE_data *data, char *filename)
//!Open the NetCDF file.
/*!
* returns  ncid: netCDF file identifier.
* @param[in] data: MLE_data struct with different MLE inputs.
* @param[in] filename: netCDF location file.
*/
{
        VERBOSE("Open the given NetCDF file... ");

        int ncid, retval;
        if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
                printf("Error: %s\n", nc_strerror(retval));
        VERBOSE("Done \n");
        return ncid;

}


int countlinesNC(int ncid, char *dim1, char *dim2)
//!Retrieve number of measurements.
/*!
* @param[in] ncid: netCDF file identifier.
* @param[in] dim1: size of dimension 1.
* @param[in] dim2: size of dimension 2.
*/
{
        int retval;
        size_t  lat_len;
        int  lat_varid;
        size_t  lon_len;
        int  lon_varid;
        //double *lats_in, *lons_in;

        //Retrieve latitiude dim length.
        if ((retval = nc_inq_dimid(ncid, dim1, &lat_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid, lat_varid, &lat_len)))
                printf("Error: %s\n", nc_strerror(retval));
        //Retrieve longitude dim length.
        if ((retval = nc_inq_dimid(ncid, dim2 , &lon_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid, lon_varid, &lon_len)))
                printf("Error: %s\n", nc_strerror(retval));

        return lat_len * lon_len;
}

void readLocsNC_1d(MLE_data *data, int ncid)
//! Read 1D locations from a given
/*! netCDF file
 * Returns location struct inside the data struct.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] ncid: netCDF file identifier.
 * */
{
        int retval;
        size_t  lat_len;
        int  lat_varid;
        size_t  lon_len;
        int  lon_varid;
        double *lats_in, *lons_in;
        int i, j, n, index;
        location *locations = &data->l1;

        VERBOSE("Read LOCS from the given NetCDF file... ");

        //Retrieve latitiude dim length.
        if ((retval = nc_inq_dimid(ncid,"latitude",&lat_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid, lat_varid, &lat_len)))
                printf("Error: %s\n", nc_strerror(retval));
        //Retrieve longitude dim length.
        if ((retval = nc_inq_dimid(ncid,"longitude", &lon_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid, lon_varid, &lon_len)))
                printf("Error: %s\n", nc_strerror(retval));

        n       = lat_len * lon_len;
        //Aloocate Memory
        lats_in         = (double *) malloc(lat_len * sizeof(double));
        lons_in         = (double *) malloc(lon_len * sizeof(double));
        locations->x     = (double *) malloc(n * sizeof(double));
        locations->y     = (double *) malloc(n * sizeof(double));

        //Get the varids of the latitude and longitude coordinate variables.
        if ((retval = nc_inq_varid(ncid, "latitude", &lat_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_varid(ncid, "longitude", &lon_varid)))
                printf("Error: %s\n", nc_strerror(retval));

        //Retrieve lats and lons values from th given file.
        if ((retval = nc_get_var_double(ncid, lat_varid, &lats_in[0])))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_get_var_double(ncid, lon_varid, &lons_in[0])))
                printf("Error: %s\n", nc_strerror(retval));

        index = 0;
        for (i = 0; i < lat_len ; i++)
                for (j = 0; j < lon_len ; j++)
                {
                        locations->x[index]      = lats_in[i];
                        locations->y[index++]    = lons_in[j];
                }

        VERBOSE("Done \n");
        //zsort_locations(n,locations);
        free(lats_in);
        free(lons_in);
}


void readLocsNC_2d(MLE_data *data, int ncid)
//! Read 2D locations from a given
/*! netCDF file
 * Returns location struct.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] ncid: netCDF file identifier.
 * */
{
        int retval;
        size_t  sn_len;
        int  lat_varid;
        size_t  we_len;
        int  lon_varid;
        //int i, j, index;
	int  n;
        location *locations = &data->l1;

        VERBOSE("Read LOCS from the given NetCDF file... ");

        //Retrieve latitiude dim length.
        if ((retval = nc_inq_dimid(ncid,"west_east",&lat_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid, lat_varid, &we_len)))
                printf("Error: %s\n", nc_strerror(retval));
        //Retrieve longitude dim length.
        if ((retval = nc_inq_dimid(ncid,"south_north", &lon_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid, lon_varid, &sn_len)))
                printf("Error: %s\n", nc_strerror(retval));

        n       = we_len * sn_len;
        //Aloocate Memory
        locations->x     = (double *) malloc(n * sizeof(double));
        locations->y     = (double *) malloc(n * sizeof(double));

        //Get the varids of the latitude and longitude coordinate variables.
        if ((retval = nc_inq_varid(ncid, "XLAT", &lat_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_varid(ncid, "XLONG", &lon_varid)))
                printf("Error: %s\n", nc_strerror(retval));

        //Retrieve lats and lons values from th given file.
        if ((retval = nc_get_var_double(ncid, lat_varid, locations->x)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_get_var_double(ncid, lon_varid, locations->y)))
                printf("Error: %s\n", nc_strerror(retval));

/*
	for(index = 0 ; index < 10; index++)
		printf("%f - %f - %f -%f\n", locations->x[index], locations->y[index]);
*/

        VERBOSE("Done \n");
        //zsort_locations(n,locations);
}

void readVarNCs(MLE_data *data, int ncid, char *varname, double *data_in , char *dim1, char *dim2)
//! Read Read observations from
/*! a  given NetCDF variable.
 * Returns location struct.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] ncid: netCDF file identifier.
 * @param[in] varname: NetCDF variable name.
 * @param[out] data_in: data vector.
 * @param[in] dim1: size of dimension 1.
 * @param[in] dim2: size of dimension 2.
 * */
{
        int varid;
        int retval;
        size_t lat_len;
        int  lat_varid;
        size_t lon_len;
        int  lon_varid;
        //int n;

        VERBOSE("Read observations from the given NetCDF variable... ");
        //Retrieve latitiude dim length.
        if ((retval = nc_inq_dimid(ncid, dim1, &lat_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid,lat_varid,&lat_len)))
                printf("Error: %s\n", nc_strerror(retval));
        //Retrieve longitude dim length.
        if ((retval = nc_inq_dimid(ncid, dim2,&lon_varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_inq_dimlen(ncid,lon_varid,&lon_len)))
                printf("Error: %s\n", nc_strerror(retval));

        //n        = lat_len * lon_len;
        size_t start[] = {0, 0, 0}; /* start at first value */
        size_t count[] = {1, lat_len, lon_len};  //wind_speed
//        size_t count[] = {1, lon_len, lat_len};  //rain, temp
        //Get the varid of the data variable, based on its name.
        if ((retval = nc_inq_varid(ncid, varname, &varid)))
                printf("Error: %s\n", nc_strerror(retval));
        if ((retval = nc_get_vara_double(ncid, varid, start, count, data_in)))
                printf("Error: %s\n", nc_strerror(retval));

        VERBOSE("DONE\n");
      
/* 
        int i;
        for(i =0;i<500;i++)
                printf("%e\n", data_in[i]);
*/
}


void closeFileNC(MLE_data *data, int ncid)
//!Close the NetCDF file.
/*!
* @param[in] ncid: netCDF file identifier.
* @param[in] data: MLE_data struct with different MLE inputs.
*/
{
        int retval;

        VERBOSE("Close the NetCDF file... ");
        if ((retval = nc_close(ncid)))
                printf("Error: %s\n", nc_strerror(retval));
        VERBOSE("Done\n");
}

