/*
 * Copyright (C) 2016, Northwestern University.
 * This file contains functions associated for cluster data structure manipulation.
 * Cluster is a double linked list structure. Kriging clustering algorithm takes advantages of this data structure for quick insertion and removal of objects.
 * To iterate through a cluster, users can see examples such as print_ functions in clusterfunctions.c
*/

#ifndef KRIG_CLUSTER_FUNCTIONS_H
#define KRIG_CLUSTER_FUNCTIONS_H

#include "cluster.h"

/*
 * Initialize a cluster.
 * Return: A new cluster readily to be used.
*/
extern Cluster* create_cluster();
/*
 * Add an object to the back of cluster.
 * cluster: cluster to be inserted to.
 * object: object being inserted.
*/
extern void add_to_cluster(Cluster* cluster,Object* object);
/*
 * Interpret a cluster as an array of objects.
 * cluster: target cluster
 * Return: An array of objects that contains all element in the cluster.
*/
extern Object** get_objects(Cluster* cluster);
/*
 * Interpret a cluster as an array of objects.
 * Extra information compared to get_objects is the size of array.
 * cluster: target cluster
 * Return: An array of objects that contains all element in the cluster.
*/
extern Objects* cluster_to_objects(Cluster* cluster);
/*
 * Clean a cluster.
 * cluster: cluster to be cleaned.
*/
extern void destroy_cluster(Cluster* cluster);
/*
 * Clean an array of clusters.
 * cluster: cluster to be destroyed.
*/
extern void destroy_clusters(Clusters* clusters);
/*
 * Print a cluster given its id to console.
 * cluster: cluster to be printed.
 * i: cluster id.
*/
extern void print_cluster(Cluster* cluster,DWORD i);
/*
 * Print all clusters to the console.
 * clusters: An array of clusters to be printed.
*/
extern void print_clusters(Clusters* clusters);
/*
 * Print objects to the console.
 * data: Pointer at the begining of the object array.
 * size: Number of objects.
*/
extern void print_objects(Object** data,DWORD size);
/*
 * Internal function used by Krig Clustering algorithm.
 * Insert an object to the middle of link list representing the cluster with O(1) complexity.
 * cluster: cluster to be inserted to.
 * previous: Previous element.
 * current: Next element.
*/
extern void insert_to_cluster(Cluster* cluster,Node* previous,Node* current);
/*
 * Internal function used by Krig Clustering algorithm.
 * Insert an object to the middle of link list representing the cluster with O(1) complexity.
 * cluster: cluster to be removed from.
 * previous: Previous element.
 * current: Next element.
*/
extern void remove_from_cluster(Cluster* cluster,Node* previous,Node* current);
/*
 * Add an object to the front of cluster.
 * cluster: cluster to be inserted to.
 * object: object being inserted.
*/
extern void add_to_cluster_front(Cluster* cluster,Object* object);
/*
 * Remove an object to the front of cluster.
 * cluster: cluster to be removed to.
*/
extern void remove_cluster_front(Cluster* cluster);
/*
 * Copy a cluster entirely.
 * cluster: The cluster to be copied.
 * Return: A new cluster tha has the same reference as the input.
*/
extern Cluster* clone_cluster(Cluster* cluster);
/*
 * Allocate memory space for training samples for variogram model.
*/
extern Samples* create_samples(DWORD size);
/*
 * Free memory space for samples.
*/
extern void destroy_samples(Samples *samples);
#endif
