/*
 * Copyright (C) 2016, Northwestern University.
 * This file contains functions for Kriging interpolations that are used by the filtering-based clustering algorithm for Kriging interpolation..
 * See also krigfunctions.h
*/
#include "matrix.h"
#include "clusterfunctions.h"
#include "krigfunctions.h"


DTYPE distance(SPATIAL_TYPE* coordinates1,SPATIAL_TYPE* coordinates2){
	DTYPE d1=coordinates1[0]-coordinates2[0];
	DTYPE d2=coordinates1[1]-coordinates2[1];
	return sqrt(d1*d1+d2*d2);
}

DTYPE exponential_variogram(DTYPE r,DTYPE C0,DTYPE C1,DTYPE C2){
	return C0+C1*(1-exp(-C2*r));
}

DTYPE power_variogram(DTYPE r,DTYPE C0,DTYPE C1,DTYPE C2){
	return C0+C1*pow(r,C2);
}

DTYPE anisotropy_power_variogram(DTYPE r,DTYPE phi,DTYPE C0,DTYPE C1,DTYPE C2,DTYPE P){
	r*=r;
	return C0+pow(pow(C1,2/P)*r*pow(cos(M_PI/4-phi),2)+pow(C2,2/P)*r*pow(cos(M_PI/4+phi),2),P/2);
}

DTYPE spherical_variogram(DTYPE r,DTYPE C0,DTYPE C1,DTYPE C2){
	return C0+2*C1*C1*(1.5*r/C2-0.5*r*r*r/(C2*C2*C2));
}

DTYPE ST_spherical_product_variogram(DTYPE r,DTYPE u,DTYPE k,DTYPE SC0,DTYPE SC1,DTYPE SC2,DTYPE TC0,DTYPE TC1,DTYPE TC2,DTYPE JC0,DTYPE JC1,DTYPE JC2){
	DTYPE intermediate;
	DTYPE vs,vt,vjoint;
	intermediate=r/SC2;
	if(intermediate<=1){
		vs=SC0+SC1*(1.5-intermediate*intermediate)*intermediate;
	}else{
		vs=0;
	}
	intermediate=u/TC2;
	if(intermediate<=1){
		vt=TC0+TC1*(1.5-intermediate*intermediate)*intermediate;
	}else{
		vt=0;
	}
	intermediate=k*u;
	intermediate=sqrt(r*r+intermediate*intermediate)/JC2;
	if(intermediate<=1){
		vjoint=JC0+JC1*(1.5-intermediate*intermediate)*intermediate;
	}else{
		vjoint=0;
	}
	return vs+vt+vjoint;
}

DTYPE ST_exp_product_variogram(DTYPE r,DTYPE u,DTYPE k,DTYPE SC0,DTYPE SC1,DTYPE SC2,DTYPE TC0,DTYPE TC1,DTYPE TC2){
	DTYPE vs=exponential_variogram(r,SC0,SC1,SC2);
	DTYPE vt=exponential_variogram(u,TC0,TC1,TC2);
	return (k+SC1)*vs+(k+TC1)*vt-k*vs*vt;
}



DTYPE compute_variogram(Object* o1,Object* o2,DTYPE* C,VARIOGRAM_TYPE variogram_type){
	DTYPE r,phi,u;
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			return exponential_variogram(distance(o1->spatial_coordinates,o2->spatial_coordinates),C[0],C[1],C[2]);
		}
		case POWER_VARIOGRAM:{
			return power_variogram(distance(o1->spatial_coordinates,o2->spatial_coordinates),C[0],C[1],C[2]);
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{
			r=distance(o1->spatial_coordinates,o2->spatial_coordinates);
			if(o1->spatial_coordinates[0]!=o2->spatial_coordinates[0]){
				phi=atan((o2->spatial_coordinates[1]-o1->spatial_coordinates[1])/(o2->spatial_coordinates[0]-o1->spatial_coordinates[0]));
			}else{
				if(o1->spatial_coordinates[0]<o2->spatial_coordinates[0]){
					phi=M_PI/2;
				}else{
					phi=-M_PI/2;
				}
			}
			return anisotropy_power_variogram(r,phi,C[0],C[1],C[2],C[3]);
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			r=distance(o1->spatial_coordinates,o2->spatial_coordinates);
			u=fabs(o1->time-o2->time);
			return ST_spherical_product_variogram(r,u,C[0],C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8],C[9]);
		}

		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
			r=distance(o1->spatial_coordinates,o2->spatial_coordinates);
			u=fabs(o1->time-o2->time);
			return ST_exp_product_variogram(r,u,C[0],C[1],C[2],C[3],C[4],C[5],C[6]);
		}
		case SPHERICAL_VARIOGRAM:{
			r=distance(o1->spatial_coordinates,o2->spatial_coordinates);
			return spherical_variogram(r,C[0],C[1],C[2]);
		}
		default:{
			return 0;
		}
	}
}
/*
 * parameters[0]=r, parameters[1]=phi, parameters[2]=u
*/
DTYPE compute_variogram_by_parameters(DTYPE *parameters,DTYPE* C,VARIOGRAM_TYPE variogram_type){
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			return exponential_variogram(parameters[0],C[0],C[1],C[2]);
		}
		case POWER_VARIOGRAM:{
			return power_variogram(parameters[0],C[0],C[1],C[2]);
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{
			return anisotropy_power_variogram(parameters[0],parameters[1],C[0],C[1],C[2],C[3]);
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			return ST_spherical_product_variogram(parameters[0],parameters[2],C[0],C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8],C[9]);
		}

		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
			return ST_exp_product_variogram(parameters[0],parameters[2],C[0],C[1],C[2],C[3],C[4],C[5],C[6]);
		}
		case SPHERICAL_VARIOGRAM:{
			return spherical_variogram(parameters[0],C[0],C[1],C[2]);
		}
		default:{
			return 0;
		}
	}
}

DTYPE DistanceSum(Cluster* cluster,Object* object,DTYPE r){
	DWORD j;
	DTYPE sum=0;
	Node* current=cluster->head;
	for(j=0;j<cluster->size;j++){
		if(current->object!=object){
			sum+=pow(distance(current->object->spatial_coordinates,object->spatial_coordinates),-r);
		}
		current=current->next;
	}
	return sum;
}

DTYPE InverseDistanceInterpolation(Cluster* cluster,Object* object,DTYPE r){
	DTYPE base=DistanceSum(cluster,object,r);
	DWORD i;
	DTYPE result=0;
	Node* current=cluster->head;
	for(i=0;i<cluster->size;i++){
		if(current->object!=object){
			result+=pow(distance(current->object->spatial_coordinates,object->spatial_coordinates),-r)*current->object->attribute;
		}
		current=current->next;
	}
	return result/base;
}

Objects* get_adjacent_objects(Cluster* cluster,Object* object,DTYPE* max_distance){
	Object** objects=Calloc(cluster->size,Object*);
	Node* front=cluster->head;
	DWORD i,index=0;
	for(i=0;i<cluster->size;i++){
		/*Either unlimited range or spatial_temporal distance within given radiuses*/
		if((max_distance[0]==DIS_UNCHECKED||distance(front->object->spatial_coordinates,object->spatial_coordinates)<max_distance[0])&&(max_distance[1]==DIS_UNCHECKED||fabs(front->object->time-object->time)<max_distance[1])){
			objects[index]=front->object;
			index++;
		}
		front=front->next;
	}
	objects=(Object**)realloc(objects,sizeof(Objects*)*index);
	Objects* result=Calloc(1,Objects);
	result->objects=objects;
	result->size=index;
	return result;
}

Matrix* krig_weights(Objects* objects,Object* object,DTYPE* C,VARIOGRAM_TYPE variogram_type){
	Object** data=objects->objects;
	Matrix* Gamma=create_matrix(objects->size+1,objects->size+1);
	Matrix* gamma=create_matrix(objects->size+1,1);
	DWORD i,j;
	for(i=0;i<objects->size;i++){
		/*Fill in the Gamma matrix, entry(i,j) is the variogram of distance between point i and point j.*/
		for(j=0;j<objects->size;j++){
			Gamma->matrix[i][j]=compute_variogram(data[i],data[j],C,variogram_type);
		}
		/*Fill in the gamma vector, each entry is the variogram of distance between the point and point i.*/
		Gamma->matrix[i][i]=0;
		gamma->matrix[i][0]=compute_variogram(object,data[i],C,variogram_type);
		//printf("data=%lf\n",data[i]->attribute);
		Gamma->matrix[i][objects->size]=1;
	}
	/*Complete filling in linear system*/
	gamma->matrix[objects->size][0]=1;
	for(j=0;j<objects->size;j++){
		Gamma->matrix[objects->size][j]=1;
	}
	Gamma->matrix[objects->size][objects->size]=0;
	//print_matrix(Gamma);
	//print_matrix(gamma);
	/*Solve the Kriging weights for the linear system*/
	Matrix* lambda=solve_linear_system(Gamma,gamma);
	/*Test code, the if condition should not be true*/
	if(lambda==NULL){
		printf("Unsolvable linear system\n");
		printf("object={%lf,%lf,%lf}\n",object->spatial_coordinates[0],object->spatial_coordinates[1],object->attribute);
		printf("first neighbor={%lf,%lf,%lf}\n",data[0]->spatial_coordinates[0],data[0]->spatial_coordinates[1],data[0]->attribute);
		printf("Gamma\n");
		print_matrix(Gamma);
		printf("gamma\n");
		print_matrix(gamma);
		return NULL;
	}
	//print_matrix(lambda);
	/*Free memory*/
	destroy_matrix(Gamma);
	destroy_matrix(gamma);
	return lambda;
}

DTYPE krig_prediction(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	DWORD i;
	Objects* objects=get_adjacent_objects(cluster,object,max_distance);
	/*If no neighbors, return the attribute without interpolation*/
	if(objects->size==0){
		return object->attribute;
	}
	Object** data=objects->objects;
	Matrix* lambda=krig_weights(objects,object,C,variogram_type);
	DTYPE result=0;
	//DTYPE sum=0;
//printf("check\n");
	/*Interpolate result using Kriging weights*/
	for(i=0;i<objects->size;i++){
		result+=lambda->matrix[i][0]*data[i]->attribute;
		//printf("data=%lf,lambda=%lf\n",data[i]->attribute,lambda->matrix[i][0]);
	//	sum+=lambda->matrix[i][0];
	}
	//print_matrix(lambda);
	//printf("sum=%lf\n",sum);

	destroy_matrix(lambda);
	Free(data);
	Free(objects);
	//printf("result=%lf\n",result);
	return result;
}
/*This function is exactly the same as the previous function except it computes Kriging variance instead of Kriging interpolation.*/
DTYPE krig_variance(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	Objects* objects=get_adjacent_objects(cluster,object,max_distance);
	if(objects->size==0){
		return 0;
	}
	Object** data=objects->objects;
	Matrix* Gamma=create_matrix(objects->size+1,objects->size+1);
	Matrix* gamma=create_matrix(objects->size+1,1);
	DWORD i,j;
	for(i=0;i<objects->size;i++){
		for(j=0;j<objects->size;j++){
			Gamma->matrix[i][j]=compute_variogram(data[i],data[j],C,variogram_type);
		}
		//Gamma->matrix[i][i]=0;
		gamma->matrix[i][0]=compute_variogram(object,data[i],C,variogram_type);
		Gamma->matrix[i][objects->size]=1;
	}
	gamma->matrix[objects->size][0]=1;
	for(j=0;j<objects->size;j++){
		Gamma->matrix[objects->size][j]=1;
	}
	Gamma->matrix[objects->size][objects->size]=0;
	Matrix* lambda=solve_linear_system(Gamma,gamma);
	destroy_matrix(Gamma);
	//printf("internal check 2\n");
	/*Compute Kriging variance from Kriging weights*/
	DTYPE result=lambda->matrix[objects->size][0];
	for(i=0;i<objects->size;i++){
		result+=(lambda->matrix[i][0])*(gamma->matrix[i][0]);
		//printf("%lf*%lf=%lf\n",lambda->matrix[i][0],gamma->matrix[i][0],(lambda->matrix[i][0])*(gamma->matrix[i][0]));
	}
	destroy_matrix(gamma);
	destroy_matrix(lambda);
	Free(data);
	Free(objects);
	return result;
}

/* 
 * This function is a combination of previous two functions.
 * Result is the normalized Kriging error.
 * Because it is more efficient to compute Kriging interpolation and Kriging variance at the same time, a separate function is provided in stead of calling these two functions separately.
 * Caching described in the implementation section of the published paper is also implemented here
*/
DTYPE krig_normalize(Cluster* cluster,Object* object,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	Objects* objects=get_adjacent_objects(cluster,object,max_distance);
	if(objects->size==0){
		return INFINITY;
		//return 0;
	}
	/*Caching*/
/*
	if(object->neighbors==objects->size){
		return object->normalized_value;
	}else{
		//printf("cache not triggered\n");
		object->neighbors=objects->size;
	}
*/
	Object** data=objects->objects;
	Matrix* Gamma=create_matrix(objects->size+1,objects->size+1);
	Matrix* gamma=create_matrix(objects->size+1,1);
	DWORD i,j;
	for(i=0;i<objects->size;i++){
		for(j=0;j<objects->size;j++){
			Gamma->matrix[i][j]=compute_variogram(data[i],data[j],C,variogram_type);
		}
		Gamma->matrix[i][i]=0;
		gamma->matrix[i][0]=compute_variogram(object,data[i],C,variogram_type);
		Gamma->matrix[i][objects->size]=1;
	}
	gamma->matrix[objects->size][0]=1;
	for(j=0;j<objects->size;j++){
		Gamma->matrix[objects->size][j]=1;
	}
	Gamma->matrix[objects->size][objects->size]=0;
	//print_matrix(gamma);
	Matrix* lambda=solve_linear_system(Gamma,gamma);
	if(lambda==NULL){
		printf("Unsolvable linear system\n");
		printf("object={%lf,%lf,%lf}\n",object->spatial_coordinates[0],object->spatial_coordinates[1],object->attribute);
		printf("first neighbor={%lf,%lf,%lf}\n",data[0]->spatial_coordinates[0],data[0]->spatial_coordinates[1],data[0]->attribute);
		printf("Gamma\n");
		print_matrix(Gamma);
		printf("gamma\n");
		print_matrix(gamma);
		return INFINITY;
	}
	destroy_matrix(Gamma);
	DTYPE result=0;
	DTYPE var=lambda->matrix[objects->size][0];
	for(i=0;i<objects->size;i++){
		result+=lambda->matrix[i][0]*data[i]->attribute;
		var+=(lambda->matrix[i][0])*(gamma->matrix[i][0]);
	}
	destroy_matrix(lambda);
	destroy_matrix(gamma);
	Free(data);
	Free(objects);
	object->normalized_value=(result-object->attribute)/sqrt(fabs(var));
	return object->normalized_value;
}
/*
 * Compute Kriging sum square error in a cluster
 * cluster: The cluster to be evaluated.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: Sum square error for leave-one-out cross validation of the cluster.
*/
DTYPE sum_krig_square_differences(Cluster* cluster,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	Node *current=cluster->head,*previous=NULL;
	DTYPE var=0,predict;
	while(current!=NULL){
		/*Leave-one-out interpolation method*/
		remove_from_cluster(cluster,previous,current);
		predict=current->object->attribute-krig_prediction(cluster,current->object,C,max_distance,variogram_type);
		var+=predict*predict;
		//DTYPE variance=krig_variance(cluster,current->object,C,max_distance,variogram_type);
		//printf("variance=%lf,normal variable=%lf\n",variance,predict);
		//printf("predicted=%lf,real=%lf\n",predict,current->object->attribute);
		insert_to_cluster(cluster,previous,current);
		previous=current;
		current=current->next;
	}
	return var;
}
/*
 * Compute the square sum of normalized Kriging error.
 * cluster: The cluster to be evaluated.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: Square sum of normalized Kriging error.
*/
DTYPE sum_krig_normalized_variance(Cluster* cluster,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	Node *current=cluster->head,*previous=NULL;
	DTYPE var=0,predict;
	while(current!=NULL){
		/*Leave-one-out interpolation method*/
		//printf("var=%lf\n",var);
		remove_from_cluster(cluster,previous,current);
		predict=krig_normalize(cluster,current->object,C,max_distance,variogram_type);
		/*
		DTYPE variance=krig_variance(cluster,current->object,C,max_distance,variogram_type);
		predict = (current->object->attribute-krig_prediction(cluster,current->object,C,max_distance,variogram_type))/sqrt(variance);
		*/
		var+=predict*predict;
		//DTYPE variance=krig_variance(cluster,current->object,C,max_distance,variogram_type);
		//printf("variance=%lf,normal variable=%lf\n",variance,predict);;
		insert_to_cluster(cluster,previous,current);
		previous=current;
		current=current->next;
	}
	return var;
}

void write_normal_squares(char* filename,Objects* objects,Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	DWORD i,j;
	Cluster *cluster,*global;
	FILE* local_copy = fopen( filename , "w" );
	fprintf(local_copy,"%s,%s\n","unclustered","clustered");
	DTYPE var,predict,benchmark;
	/*For each clusters*/
	for(i=0;i<clusters->size;i++){
		if(clusters->clusters[i]->size<2){
			continue;
		}
		cluster=clusters->clusters[i];
		Node *current=cluster->head,*previous=NULL;
		/*For each element in the cluster*/
		while(current!=NULL){
			/*Write normalized Kriging error for both Kriging with clustering and Kriging without clustering*/
			remove_from_cluster(cluster,previous,current);
			predict=krig_prediction(cluster,current->object,C,max_distance,variogram_type);
			var=fabs(predict-current->object->attribute);
			global=create_cluster();
			for(j=0;j<objects->size;j++){
				if(objects->objects[j]!=current->object){
					add_to_cluster(global,objects->objects[j]);
				}
			}
			benchmark=krig_prediction(global,current->object,C,max_distance,variogram_type);
			destroy_cluster(global);
			benchmark=fabs(benchmark-current->object->attribute);
			if(benchmark!=INFINITY&&var!=INFINITY){
				fprintf(local_copy,"%lf,%lf\n",benchmark,var);
			}
			insert_to_cluster(cluster,previous,current);
			previous=current;
			current=current->next;
		}
	}
}
/*
 * Compute the chi-square statistics for an array of clusters.
 * clusters: The array of clusters.
 * C: Parameters for Kriging interpolation.
 *    Exponential variogram: C has size 3. Variogram is calculated as C[0]+C[1]*(1-exp(-C[2]*r)).
 *    Power variogram: C has size 3. Variogram is calculated as C[0]+C[1]*pow(r,C[2])
 *    Anisotropy_power_variogram: C has size 4. Variogram is calculated as C[0]+pow(pow(C[1],2/C[3])*r*pow(cos(M_PI/4-phi),2)+pow(C[2],2/C[3])*r*pow(cos(M_PI/4+phi),2),C[3]/2)
 *    Spatio-temporal spherical product variogram: C has size 10. Variogram is calculated as vs+vt+vjoint. vs=C[1]+C[2]*(1.5-(r/C[3])^2)*r/C[3], vt=C[4]+C[5]*(1.5-(r/C[6])^2)*t/C[6], and vjoint=C[7]+C[8]*(1.5-(sqrt(s^2+(C[0]*t)^2)/C[9])^2)*sqrt(s^2+(C[0]t)
 *    Spatio-temporal exponential product variogram: C has size 7. (C[0]+C[2])*vs+(C[0]+C[5])*vt+C[0]vs*vt. vs=C[1]+C[2]*pow(r,C[3]), vt=C[4]+C[5]*pow(t,C[6])
 * max_distance: Array of size 2. The first element is the maximum spatial ball range for local Kriging. The second element is maximum temporal ball range for local Kriging. If local Kriging is not required, these values should be set to constant DIS_UNCHECKED.
 * variogram_type: The type of variogram. Constant values in "clustertype.h"
 * Return: Chi-square statistics for the array of clusters cluster (sum of square of all nomalized Kriging errors)
*/
DTYPE chi_square_coefficient(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	DTYPE x=0;
	DTYPE temp;
	DWORD i;
	for(i=0;i<clusters->size;i++){
		//printf("cluster=%lld,x=%lf\n",i,x);
		if(clusters->clusters[i]->size>1){
			temp=sum_krig_normalized_variance(clusters->clusters[i],C,max_distance,variogram_type);
			x+=temp;
		}
	}
	//printf("\n");
	return x;
}
/*
 * Compute the variance of physical attributes for all points in a cluster.
 * cluster: The cluster to be evaluated.
 * Return: The mean of physical attribute for all points in the cluster.
*/
DTYPE cluster_mean(Cluster* cluster){
	DTYPE result=0;
	DWORD i;
	Node* front=cluster->head;
	for(i=0;i<cluster->size;i++){
		result+=front->object->attribute;
		front=front->next;
	}
	return result/cluster->size;
}
/*
 * Compute the variance of physical attributes for all points in a cluster.
 * cluster: The cluster to be evaluated.
 * Return: The variance of physical attribute for all points in the cluster.
*/
DTYPE cluster_variance(Cluster* cluster){
	DTYPE mean=cluster_mean(cluster);
	DTYPE result=0;
	DTYPE var;
	DWORD i;
	Node* front=cluster->head;
	for(i=0;i<cluster->size;i++){
		var=(front->object->attribute-mean);
		var*=var;
		result+=var;
		front=front->next;
	}
	return result/(cluster->size);
}
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
DTYPE square_errors(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	DTYPE x=0;
	DWORD i;
	for(i=0;i<clusters->size;i++){
		//printf("cluster=%lld,x=%lf\n",i,x);
		if(clusters->clusters[i]->size>1){
			x+=sum_krig_square_differences(clusters->clusters[i],C,max_distance,variogram_type);
		}
	}
	return x;
}

DTYPE NMSE_error(Clusters* clusters,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	DTYPE x=0;
	DWORD i;
	for(i=0;i<clusters->size;i++){
		//printf("cluster=%lld,x=%lf\n",i,x);
		if(clusters->clusters[i]->size>10){
			//printf("variance*n=%lf\n",cluster_variance(clusters->clusters[i])*clusters->clusters[i]->size);
			//printf("SE=%lf\n",sum_krig_square_differences(clusters->clusters[i],C,max_distance,variogram_type));
			x+=sum_krig_square_differences(clusters->clusters[i],C,max_distance,variogram_type)/(cluster_variance(clusters->clusters[i])*clusters->clusters[i]->size);
			//printf("term=%lf\n",x);
		}
	}
	return x;
}
/*Used for sorting points based on their normalized Kriging error*/
int object_cmp(const void* o1,const void* o2){
	DTYPE d1=(*((Object**)o1))->normalized_value;
	DTYPE d2=(*((Object**)o2))->normalized_value;
	return d1<d2?1:-1;
}

BOOLEAN krig_consistency(Cluster* cluster,DTYPE bound,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	Object** clone=get_objects(cluster);
	/*Heuristic implementation discussed in the imlementation section of the published paper.*/
	qsort(clone,cluster->size,sizeof(Object*),object_cmp);
	DWORD i;
	Cluster* copy=create_cluster();
	for(i=0;i<cluster->size;i++){
		add_to_cluster(copy,clone[i]);
	}
	DTYPE predict;
	Node *current=copy->head,*previous=NULL;
	while(copy->size>1&&current!=NULL){
		//printf("check\n");
		remove_from_cluster(copy,previous,current);
		predict=krig_normalize(copy,current->object,C,max_distance,variogram_type);
		insert_to_cluster(copy,previous,current);
		if(fabs(predict)>bound){
			return FALSE;
		}
		previous=current;
		current=current->next;
	}
	destroy_cluster(copy);
	Free(clone);
/*
	Node *current=cluster->head,*previous=NULL;
	while(cluster->size>1&&current!=NULL){
		//printf("check\n");
		remove_from_cluster(cluster,previous,current);
		predict=krig_normalize(cluster,current->object,C,max_distance,variogram_type);
		insert_to_cluster(cluster,previous,current);
		if(fabs(predict)>bound){
			return FALSE;
		}
		previous=current;
		current=current->next;
	}
*/
	return TRUE;
}

DTYPE sum_krig_variance(Cluster* cluster,DTYPE* C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	Node *current=cluster->head,*previous=NULL;
	DTYPE var=0;
	while(current!=NULL){
		/*Leave-one-out interpolation method*/
		remove_from_cluster(cluster,previous,current);
		var+=krig_variance(cluster,current->object,C,max_distance,variogram_type);
		insert_to_cluster(cluster,previous,current);
		previous=current;
		current=current->next;
	}
	return var;
}
