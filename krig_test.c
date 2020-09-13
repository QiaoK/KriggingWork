/*
*  Copyright (C) 2016, Northwestern University.
*/
#include "cluster.h"
#include "clusterfunctions.h"
#include "krigfunctions.h"
/*
 * This file is an example for how to use Kriging prediction functions in this software package.
 * This function is a testing funtion for Kriging interpolation.
 * The first example used is in Sherman's book chapter 2.2.
 * The second example is a robust test for repeated spatial coordinates.
*/


int main(void){
	//Fill up all objects
	Object** objects=Calloc(3,Object*);
	objects[0]=Calloc(1,Object);
	objects[1]=Calloc(1,Object);
	objects[2]=Calloc(1,Object);
	objects[0]->spatial_coordinates=Calloc(2,DTYPE);
	objects[1]->spatial_coordinates=Calloc(2,DTYPE);
	objects[2]->spatial_coordinates=Calloc(2,DTYPE);
	objects[1]->attribute=2.412;
	objects[2]->attribute=4.525;
	objects[1]->spatial_coordinates[0]=0;
	objects[1]->spatial_coordinates[1]=0;
	objects[2]->spatial_coordinates[0]=3+47.0/6;
	objects[2]->spatial_coordinates[1]=sqrt(169-47.0*47/36);
	objects[0]->spatial_coordinates[0]=3;
	objects[0]->spatial_coordinates[1]=0;
	printf("d(s1,s2)=%lf\n",distance(objects[1]->spatial_coordinates,objects[2]->spatial_coordinates));
	//Test for cluster function
	Cluster* cluster=create_cluster();
	add_to_cluster(cluster,objects[1]);
	add_to_cluster(cluster,objects[2]);
	print_cluster(cluster,0);
	//Construct variogram function
	DTYPE* C=Calloc(2,DTYPE);
	C[0]=.5;
	C[1]=.2;
	DTYPE* T=Calloc(2,DTYPE);
	T[0]=DIS_UNCHECKED;
	T[1]=DIS_UNCHECKED;
	//Interpolate result
	DTYPE prediction=krig_prediction(cluster,objects[0],C,T,EXPONENTIAL_VARIOGRAM);
	//Compute variance
	DTYPE variance=krig_variance(cluster,objects[0],C,T,EXPONENTIAL_VARIOGRAM);
	//Print result
	printf("prediction=%lf,variance=%lf\n",prediction,variance);
	
	//Test for duplication
	objects[1]->spatial_coordinates[0]=0;
	objects[1]->spatial_coordinates[1]=0;
	objects[2]->spatial_coordinates[0]=0;
	objects[2]->spatial_coordinates[1]=0;
	objects[0]->spatial_coordinates[0]=3;
	objects[0]->spatial_coordinates[1]=0;
	prediction=krig_prediction(cluster,objects[0],C,T,EXPONENTIAL_VARIOGRAM);
	variance=krig_variance(cluster,objects[0],C,T,EXPONENTIAL_VARIOGRAM);
	printf("prediction=%lf,variance=%lf\n",prediction,variance);

	Free(objects[0]->spatial_coordinates);
	Free(objects[1]->spatial_coordinates);
	Free(objects[2]->spatial_coordinates);
	Free(objects);
	return 0;
}
