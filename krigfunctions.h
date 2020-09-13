/*
 * Copyright (C) 2016, Northwestern University.
 * User visible interfaces.
 * All functions for Kriging clustering and Kriging interpolation methods are listed below.
 * See krig_functions.c and krig_cluster.c
*/

#ifndef KRIG_FUNCTIONS_H
#define KRIG_FUNCTIONS_H

#include "cluster.h"
#include "matrix.h"

/*
 * Compute distance between two spatial coordinates. 2D data is assumed.
 * Modify the function in krig_functions.c for high dimensional data usage.
 * coordiantes1: Spatial coordinate of first point.
 * coordinates2: Spatial coordinate of second point.
*/
extern DTYPE distance(SPATIAL_TYPE* coordinates1,SPATIAL_TYPE* coordinates2);
/*
 * Compute distance between an object and every element in a cluster. Sum the distances to the power of r.
 * Exclusively used for inverse distance interpolation, not Kriging interpolation.
 * cluster: The cluster that are used to compute sum distance to power r towards the object.
 * object: The object that are used to compute sum distance to power r towards the cluster.
 * r: power. e.g for sum squares, r=2. For linear sum, r=1.
 * Return: Sum the distances to the power of r.
*/
extern DTYPE DistanceSum(Cluster* cluster,Object* object,DTYPE r);
/*
 * Interpolate the physical attribute of an object using all objects in a cluster by inverse distance interpolation method.
 * cluster: The cluster which contains objects that are used to interpolate the attribute of the object.
 * object: The object which its physical attribute is interpolated.
 * r: The power used for inversed distance interpolation.
*/
extern DTYPE InverseDistanceInterpolation(Cluster* cluster,Object* object,DTYPE r);
/*
 * Interpolate the physical attribute of an object using all objects in a cluster by clustering-based Kriging interpolation algorithm.
 * cluster: The cluster which contains objects that are used to interpolate the attribute of the object.
 * object: The object which its physical attribute is interpolated.
 * r: The power used for inversed distance interpolation.
*/
extern DTYPE krig_prediction(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Filtering-based Kriging clustering function.
 * Proposed in "A Filtering-based Clustering Algorithm for Improving Spatio-temporal Kriging Interpolation Accuracy", CIKM 2016
 * data: An array of objects that will be clustered.
 * size: The number of objects.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: An array of clusters satisfying consistency and maximality constraints discussed in the paper.
*/
extern Clusters* krig_clustering(Object** data,DWORD size,DTYPE bound,DTYPE *C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Evaluate a set of clusters by computing the chi-square coefficient.
 * This measurement was proposed in "A Filtering-based Clustering Algorithm for Improving Spatio-temporal Kriging Interpolation Accuracy", CIKM 2016
 * Result is the square sum of normalized clustering-based Kriging interpolation error.
 * clusters: The set of clusters being evaluated.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: The chi-square error.
*/
extern DTYPE chi_square_coefficient(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Compute the normalized Kriging error for an object using clustering-based Kriging interpolation by (prediction-actual)/variance.
 * cluster: The cluster used for clustering-based Kriging interpolation.
 * object: The object which its physical attribute will be evaluates.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: The normalized Kriging error.
*/
extern DTYPE krig_normalize(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Check if a cluster is consistent in terms of normalized clustering-based Kriging interpolation error.
 * Definition of consistency was proposed in "A Filtering-based Clustering Algorithm for Improving Spatio-temporal Kriging Interpolation Accuracy", CIKM 2016
 * cluster: The cluster to be evaluated.
 * bound: The threshold for consistency.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: If the cluster is consistent.
*/
extern BOOLEAN krig_consistency(Cluster* cluster,DTYPE bound,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Compute the clustering-based Kriging variance for an object's physical attribute.
 * cluster: The cluster used for clustering-based Kriging.
 * object: The object that has the physical attribute for variance computation.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])Samples* variogram_sampling(Objects *objects,DTYPE bound,DTYPE angle_bound,DTYPE step_size,DWORD steps,DTYPE *C,SMOOTHING_TYPE smoothing_type)
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: Clustering-based Kriging variance
*/
extern DTYPE krig_variance(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Evaluate a set of clusters by computing the normalized mean square error.
 * clusters: The set of clusters being evaluated.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: The NMSE error.
*/
extern DTYPE NMSE_error(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * Get spatio-temporal neighbors of an object in a cluster.
 * cluster: The cluster that contains potential neighbors.
 * object: The object which its neighbors will be found.
 * max_distance: Array of size 2. The first element is the maximum spatial radius for neighbors. The second element is maximum temporal radius for neighbors.
*/
extern Objects* get_adjacent_objects(Cluster* cluster,Object* object,DTYPE* max_distance);
/*
 * Compute the sum square errors for an array of clusters.
 * clusters: The array of clusters.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: Sum square errors for an array of clusters. (Sum of individual cluster square error.)
*/
extern DTYPE square_errors(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * This function write the result of absolute Kriging error for Kriging with clustering and Kriging without clustering to a local CSV file.
 * filename: The file to write to.
 * objects: The array of all objects. This input is required for constructing Kriging without clustering scenario.
 * clusters: The clusters for clustering-based Kriging interpolation.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
*/
extern void write_normal_squares(char* filename,Objects* objects,Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
/*
 * This function computes estimation of variograms based on input parameters.
 * parameters: 0 is spatial lag h, 1 is spatial angle phi, 2 is temporal lag u.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
*/
extern DTYPE sum_krig_square_differences(Cluster* cluster,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
extern DTYPE sum_krig_normalized_variance(Cluster* cluster,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
extern DTYPE sum_krig_variance(Cluster* cluster,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);
extern Matrix* krig_weights(Objects* objects,Object* object,DTYPE* C,VARIOGRAM_TYPE variogram_type);

//Variogram related functions.
extern DWORD variogram_model_length(VARIOGRAM_TYPE variogram_type);
extern DTYPE compute_variogram_by_parameters(DTYPE *parameters,DTYPE* C,VARIOGRAM_TYPE variogram_type);
extern Samples* variogram_sampling(Objects *objects,DTYPE bound,DTYPE angle_bound,DTYPE step_size,DWORD steps,DWORD angle_steps,DTYPE *C,SMOOTHING_TYPE smoothing_type);
extern DTYPE evaluate_model(Samples* samples, DTYPE* C, VARIOGRAM_TYPE variogram_type);
extern DTYPE* variogram_PSO(Objects* objects,Samples* samples,DWORD epochs, DWORD n_particles, VARIOGRAM_TYPE variogram_type, DTYPE c1,DTYPE c2,DTYPE alpha);
extern void variogram_WLS(Samples* samples,DTYPE* C,DWORD epochs, VARIOGRAM_TYPE variogram_type,DTYPE learning_rate);
extern void variogram_WLS_line_search(Samples* samples,DTYPE *C, DWORD epochs, VARIOGRAM_TYPE variogram_type);
extern void randomize_variogram(Samples* samples, DTYPE* C, VARIOGRAM_TYPE variogram_type);
extern void print_samples(Samples* samples);
#endif
