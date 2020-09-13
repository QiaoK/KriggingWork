/*
 * Copyright (C) 2016, Northwestern University.
 * This file contains essential definitions for data structures that are used in this pacakge.
 * See also clustertype.h for internal type and clusterfunctions.h for operations on data structure.
*/

#ifndef KRIG_CLUSTER_H
#define KRIG_CLUSTER_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include "clustertype.h"

/*
 * Type definition for basic object structure.
 * time : The time stamp that the object is located.
 * attribute: The physical attribute of this object.
 * spatial_coordinates: The spatial coordinate that the object is at.
 * normalized_valkue: Caching value for current normalized Kriging error. Only used by Kriging clustering algorithm at low level.
 * neighbors: Number of spatio-temporal neighbors around this object. Only used by Kriging clustering algorithm at low level.
*/
typedef struct{
	DTYPE time;
	DTYPE attribute;
	SPATIAL_TYPE* spatial_coordinates;
	DTYPE normalized_value;
	DWORD neighbors;
} Object;

/*
 * Link list for objects;
*/
typedef struct Node{
	Object* object;
	struct Node *next;
} Node;

/*
 * An array of objects.
*/
typedef struct{
	Object** objects;
	DWORD size;
} Objects;

/*
 * A data structure for cluster that contains a double linked list of objects;
 * Sum of physical attributes are computed at run time.
*/

typedef struct{
	Node* head;
	DWORD size;
	DTYPE sum;
	Node* tail;
} Cluster;

/*
 * An list of objects representing neighbors of an object. Exactly the same as Objects.
*/
typedef struct{
	Object** neighbors;
	DWORD size;
} Neighbors;

/*
 * An array of clusters
*/
typedef struct{
	Cluster** clusters;
	DWORD size;
} Clusters;
/*
 * Training sample for variogram
*/
typedef struct{
	DTYPE* x;
	DTYPE* phi;
	DTYPE* y;
	DWORD* N;
	DWORD size;
} Samples;
/*
 * Particle used in PSO algorithm
*/
typedef struct{
	DTYPE* C;
	DTYPE* velocity;
	DTYPE* local_best;
	DTYPE mse;
} Particle;

#endif
