/*
*  Copyright (C) 2016, Northwestern University.
*/
#include "cluster.h"
#include "krigfunctions.h"
#include "random.h"
#include "clusterfunctions.h"
#include "datafunctions.h"
/*
 * This file explains how to use Kriging clustering algorithm for SOCR data
*/
void print_data(Objects* data){
	DWORD i;
	for(i=0;i<data->size;i++){
		printf("x=%lf,y=%lf,attribute=%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->attribute);
	}
}
/*
 * Read in K-means clustering result in format produced by
 * http://users.eecs.northwestern.edu/~wkliao/Kmeans/
 * Number of clusters is assumed to be 6.
*/
Clusters* read_kmeans(Objects* data){
	Clusters* clusters=Calloc(1,Clusters);
	clusters->clusters=Calloc(6,Cluster*);
	clusters->size=6;
	//Double stream for special character capturing
	FILE* stream1=fopen("SOCR.txt.membership","r");
	FILE* stream2=fopen("SOCR.txt.membership","r");
	DWORD i;
	for(i=0;i<6;i++){
		clusters->clusters[i]=create_cluster();
	}
	DWORD length=0;
	char* line,*temp;
	char c;
	i=0;
	while((c=fgetc(stream1))!=EOF){
		length++;
		if(c=='\n'&&length>=1){
			line=Calloc(length+1,char);
			fgets(line,length,stream2);
			//printf("%s\n",line);
			fgetc(stream2);
			line[length]='\0';
			temp=line;
			while(*temp!=' '){
				temp++;
			}
			temp++;
			//printf("processing %d\n",atoi(temp));
			add_to_cluster(clusters->clusters[atoi(temp)],data->objects[i]);
			i++;
			length=0;
			Free(line);
		}
	}
	fclose(stream1);
	fclose(stream2);
	return clusters;
}
/*
 * Test for SOCR data. NMSE, chi-square and number of clusters will be printed. Timing information will also be printed. Each cluster will be written into an individual file that shows the objects in the cluster. For example, SOCR_0.txt contains all objects in cluster with id 0. Comparison of clustered statistics and unclustered statistics will be written into SOCR_compare.txt
*/
int main(void){
	//Read SOCR data
	Objects* data=read_csv("SOCR_061708_NC_Data_Aquifer.csv",SPATIAL_DATA,FALSE);
	DWORD i;
	//Rescaling data
	for(i=0;i<data->size;i++){
		data->objects[i]->attribute*=1000;
	}
	//shuffle_objects(data->objects,data->size);
	//print_data(data);
	//Construct variogram
	DTYPE* C=Calloc(4,DTYPE);
	C[0]=14000;
	C[1]=38;
	C[2]=15;
	C[3]=1.99;
	DTYPE* distances=Calloc(2,DTYPE);
	distances[0]=DIS_UNCHECKED;
	distances[1]=DIS_UNCHECKED;
	//Kriging clustering for data.
	Clusters* clusters=krig_clustering(data->objects,data->size,1,C,distances,ANISOTROPHY_POWER_VARIOGRAM);
	print_clusters(clusters);
	//Write clusters to local file.
	write_clusters(clusters,"SOCR_",TRUE);
	write_normal_squares("SOCR_compare.txt",data,clusters,C,distances,ANISOTROPHY_POWER_VARIOGRAM);
	//Print Chi-square result
	DTYPE test_statistics=chi_square_coefficient(clusters,C,distances,ANISOTROPHY_POWER_VARIOGRAM);
	printf("chi square test statistics=%lf\n",test_statistics);
	//Print NMSE
	printf("NMSE error=%lf\n",NMSE_error(clusters,C,distances,ANISOTROPHY_POWER_VARIOGRAM));
	//Print results produced by K-means clustering.
	Clusters* kmeans=read_kmeans(data);
	write_normal_squares("SOCR_kmeans_compare.txt",data,kmeans,C,distances,ANISOTROPHY_POWER_VARIOGRAM);
	printf("NMSE error=%lf\n",NMSE_error(kmeans,C,distances,ANISOTROPHY_POWER_VARIOGRAM));
	test_statistics=chi_square_coefficient(kmeans,C,distances,ANISOTROPHY_POWER_VARIOGRAM);
	printf("chi square test statistics=%lf\n",test_statistics);
	//print_clusters(kmeans);
	Free(C);
	Free(distances);
	destroy_clusters(clusters);
	destroy_clusters(kmeans);
	return 0;
}
