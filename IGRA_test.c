/*
*  Copyright (C) 2016, Northwestern University.
*/
#include "krigfunctions.h"
#include "random.h"
#include "clustertype.h"
#include "cluster.h"
#include "clusterfunctions.h"
#include "datafunctions.h"
/*
 * Example for how to use Kriging clustering algorithm on IGRA dataset.
 * 60 days data are read in. Chi-square values for Kriging with clustering and Kriging without clustering are computed and write to csv files for each time stamps. For example IGRA_compare_0.txt is comparison for time stamp 0.
*/
int main(void){
	//Read IGRA data
	Objects* data=read_IGRA("US_IGRA",60,221,TRUE);
	DWORD i,j,index=0;
	char file_name[100];
	//printf("filtered size=%lld\n",index);
	//Construct variogram
	DTYPE* C=Calloc(3,DTYPE);
	C[0]=54.84102;
	C[1]=153.71639;
	C[2]=0.05462576;
	DTYPE* distances=Calloc(2,DTYPE);
	distances[0]=30;
	distances[1]=DIS_UNCHECKED;
	DTYPE increment=-88;
	for(j=0;j<60;j++){
		//Data for day i are read
		Object** filtered=Calloc(data->size,Object*);
		index=0;
		if(j%2==0){
			increment+=88;
		}else{
			increment+=12;
		}
		for(i=0;i<data->size;i++){
			//printf("x=%lf,y=%lf,t=%lf,time=%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->attribute,data->objects[i]->time);
			if(data->objects[i]->time==2014010100+increment&&data->objects[i]->attribute!=-9999){
				filtered[index++]=data->objects[i];
			}
		}
		//printf("index=%lld\n",index);
		filtered=(Object**)realloc(filtered,sizeof(Object*)*index);
		//shuffle_objects(filtered,index);
		//Kriging without clustering (unlimited error boundary, everything will be in the same cluster)
		Clusters* clusters=krig_clustering(filtered,index,99999,C,distances,EXPONENTIAL_VARIOGRAM);
		//Compute chi-sqaure value for Kriging without clustering.
		DTYPE benchmark=chi_square_coefficient(clusters,C,distances,EXPONENTIAL_VARIOGRAM);
		destroy_clusters(clusters);
		//Kriging with clustering, threshold is 0.6.
		clusters=krig_clustering(filtered,index,0.6,C,distances,EXPONENTIAL_VARIOGRAM);
		//print_clusters(clusters);
		//Compute chi-sqaure value for Kriging with clustering
		DTYPE test_statistics=chi_square_coefficient(clusters,C,distances,EXPONENTIAL_VARIOGRAM);
		//write out day i's comparison result
		Objects* objects=Calloc(1,Objects);
		objects->objects=filtered;
		objects->size=index;
		sprintf(file_name, "IGRA_compare_%lld.txt",j);
		write_normal_squares(file_name,objects,clusters,C,distances,EXPONENTIAL_VARIOGRAM);
		printf("%lld,%lld,%lf,%lf,%lf\n",index,clusters->size,benchmark,test_statistics,NMSE_error(clusters,C,distances,EXPONENTIAL_VARIOGRAM));
		Free(filtered);
		Free(objects);
	}
	Free(C);
	Free(distances);
	return 0;
}
