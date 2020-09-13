This software implements filtering-based clusterings algorithm for improving Kriging interpolation accuracy.
Main parts: Matrix operations, Kriging interpolations, Clustering functions, and filtering-based clustering.
Input and output files are all text files.

# System requirement

gcc compiler that supports C89.

# Compile

Use make file to compile code and test cases with command "make".


# Data types

Basic types and constants are defined in clustertype.h. For example, DTYPE means double (floating point). DWORD means long long int (integer).
Structure types are defined in cluster.h.


# Functions

These functions are interfaces that users may be interested in. Include krigfunctions.h in order to use them.

/*
 * Interpolate the physical attribute of an object using all objects in a cluster by clustering-based Kriging interpolation algorithm.
 * cluster: The cluster which contains objects that are used to interpolate the attribute of the object.
 * object: The object which its physical attribute is interpolated.
 * r: The power used for inversed distance interpolation.
*/
DTYPE krig_prediction(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);

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
Clusters* krig_clustering(Object** data,DWORD size,DTYPE bound,DTYPE *C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);

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
DTYPE krig_normalize(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);

/*
 * Compute the clustering-based Kriging variance for an object's physical attribute.
 * cluster: The cluster used for clustering-based Kriging.
 * object: The object that has the physical attribute for variance computation.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: Clustering-based Kriging variance
*/
DTYPE krig_variance(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);

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
DTYPE NMSE_error(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type);


# Executables

All executables have built in inputs without arguments.

- ./IGRA_test	
 60 days data are read in. Chi-square values for Kriging with clustering and Kriging without clustering are computed and write to csv files for each time stamps. For example IGRA_compare_0.txt is comparison for time stamp 0.

- ./SOCR_test
 Test for SOCR data. NMSE, chi-square and number of clusters will be printed. Timing information will also be printed. Each cluster will be written into an individual file that shows the objects in the cluster. For example, SOCR_0.txt contains all objects in cluster with id 0. Comparison of clustered statistics and unclustered statistics will be written into SOCR_compare.txt

- ./matrix_test
 Unit test for matrix operations. Should pass all test cases without any alerts.

# Contact

Contact me at qiao.kang@eecs.northwestern.edu if you have questions.
