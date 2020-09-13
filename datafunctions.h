/*
 * Copyright (C) 2016, Northwestern University.
 * This file contains data I/O functions that are necessary for processing test datasets.
 * Test datasets include SOCR, IGRA, and OZONE tram. Input datafiles are assumed to be in text format.
 * See also SOCR_test, IGRA_test, and pollution_test for example usage.
*/
#ifndef KRIG_DATA_FUNCTION_H
#define KRIG_DATA_FUNCTION_H

#include "cluster.h"

//Spatial data identifier
#define SPATIAL_DATA 0x1415
//Spatio-temporal data identifier
#define SPATIAL_TEMPORAL_DATA 0x1963


/*
 * Read raw IGRA data from a folder and convert the data into objects.
 * dir_name: The directory that contains all IGRA data.
 * time_elapse: How many distinct time stamps a user want to read from these data.
 * station_size: Number of spatial coordinates, which is the number of folders that are contained in the folder.
 * local_copy_flag: If you want to make a csv copy of IGRA data for future use.
 * Return: An array of objects for IGRA data.
*/
extern Objects* read_IGRA(char* dir_name,DWORD time_elapse,DWORD station_size,BOOLEAN local_copy_flag);
/*
 * Read a csv file that either contain spatial data or spatio-temporal data.
 * For case of spatial data, columns must be (x, y, attribute)
 * For case of spatio-temporal data, columns must be (x, y, attribute, time stamp)
 * filename: The file to be read.
 * type: SPATIAL_DATA or SPATIAL_TEMPORAL_DATA are accepted for two types of data.
*/
extern Objects* read_csv(char* filename,WORD type,BOOLEAN header);
/*
 * Write spatio-temporal objects to local csv file.
 * Columns are (x, y, attribute, time stamp)
 * data: Array of objects that will be exported.
 * filename: The file to be written.
 * header: If headers for columns should be added.
*/
extern void write_spatial_temporal_data(Objects* data,char* filename,BOOLEAN header);
/*
 * Write an array of clusters to local file system.
 * Each cluster will be represented by a single file with in form 'filename'_'clusterid'.txt
 * Objects are represented in column form (x,y,time,attribute).
 * clusters: The array of clusters that will be exported.
 * header: If headers for columns should be added.
*/
extern void write_clusters(Clusters* clusters,char* filename,BOOLEAN header);

extern void write_variogram_result(char* filename, BOOLEAN header, Samples* samples, DTYPE* C, VARIOGRAM_TYPE variogram_type);
#endif
