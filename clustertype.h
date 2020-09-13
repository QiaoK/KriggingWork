/*
 * Copyright (C) 2016, Northwestern University.
 * Internal type used in this package. Changing basic type definition will be valid for the entire software package.
 * e.g double to float will automatically reset the precision of all floating point operations.
*/
#ifndef KRIG_CLUSTER_TYPE_H
#define KRIG_CLUSTER_TYPE_H

//Basic type for integer
#define DWORD long long
//Basic type for small range integer type
#define WORD int
//Type for floating point, default precision is double
#define DTYPE double
//Type for spatial coordinate
#define SPATIAL_TYPE double
//Type for attribute
#define NON_SPATIAL_TYPE double
//Type for time stamp
#define TEMPORAL_TYPE double
//Type for kernel smoother
#define SMOOTHING_TYPE WORD
//Gaussian smoothing
#define GAUSSIAN_SMOOTHING 0x1211
//Uniform smoothing
#define UNIFORM_SMOOTHING 0x1212
//Type for variogram
#define VARIOGRAM_TYPE WORD
//Exponential variogram constant value
#define EXPONENTIAL_VARIOGRAM 0x2513
//Power variogram constant value
#define POWER_VARIOGRAM 0x4353
//Anisotrophy power variogram constant value
#define ANISOTROPHY_POWER_VARIOGRAM 0x3235
//Spherical variogram constant value
#define SPHERICAL_VARIOGRAM 0x3236
//Spatio-temporal spherical product variogram constant value
#define ST_SPHERICAL_PRODUCT_VARIOGRAM 0x7235
//Spatio-temporal exponential product variogram constant value
#define ST_EXPONENTIAL_PRODUCT_VARIOGRAM 0x7236
//Type for temporal distance
#define TEMPORAL_DISTANCE 0x725
//Type for spatial distance
#define SPATIAL_DISTANCE 0x772
//Boolean values
#define TRUE 1
#define FALSE 0
#define BOOLEAN char
/*
 * Memory allocator. stdlib.h used. Replace with your own memory allocation function.
 * a: number of elements.
 * b: type of element.
*/
#define Calloc(a,b) ((b*)malloc(sizeof(b)*a))
#define Free free
#define DIS_UNCHECKED -1
#endif
