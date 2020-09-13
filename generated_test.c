#include "krigfunctions.h"
#include "random.h"
#include "clusterfunctions.h"
#include "datafunctions.h"
/*
 * Kriging clustering example for random data.
*/


void print_data(Objects* data){
	DWORD i;
	for(i=0;i<data->size;i++){
		printf("x=%lf,y=%lf,attribute=%lf,time=%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->attribute,data->objects[i]->time);
	}
}

int main(void){
/*
	struct timeval start;
	gettimeofday(&start,NULL);
	init_genrand(start.tv_usec);
	DWORD cluster_size=10,max_size=500,min_points=10;
	DTYPE max_distance=0.02,max_distance2=2.0,x_min=0,x_max=1,y_min=0,y_max=1,t_min=30,t_max=50,alpha=1.3;
	Stack* stack=random_data(cluster_size,max_size,min_points,max_distance,max_distance2,x_min,x_max,y_min,y_max,t_max,t_min,alpha);
	Object** objects=Calloc(max_size,Object*);
	DWORD i=0;
	while(!is_empty(stack)){
		objects[i]=pop(stack);
		objects[i]->classification=UNCLASSIFIED;
		i++;
	}
	DTYPE* C=Calloc(2,DTYPE);
	C[0]=.2;
	C[1]=.5;
	DTYPE* distances=Calloc(2,DTYPE);
	distances[0]=max_distance2;
	distances[1]=DIS_UNCHECKED;
	Objects* data=Calloc(1,Objects);
	data->size=max_size;
	data->objects=objects;
	write_spatial_temporal_data(data,"random.csv",TRUE);
*/
	Objects* data=read_csv("random.csv",SPATIAL_TEMPORAL_DATA,TRUE);
	print_data(data);
	Object** objects=data->objects;
	DWORD cluster_size=10,max_size=500,min_points=10;
	DTYPE max_distance=0.02,max_distance2=.3,x_min=0,x_max=1,y_min=0,y_max=1,t_min=30,t_max=50,alpha=1.3;
	DTYPE* C=Calloc(2,DTYPE);
	C[0]=0.4;
	C[1]=.4;
	DTYPE* distances=Calloc(2,DTYPE);
	distances[0]=max_distance2;
	distances[1]=DIS_UNCHECKED;
	Clusters* clusters=krig_clustering(objects,max_size,3,C,distances,EXPONENTIAL_VARIOGRAM);
	printf("# of clusters=%lld\n",clusters->size);
	print_clusters(clusters);
	DTYPE test_statistics=chi_square_coefficient(clusters,C,distances,EXPONENTIAL_VARIOGRAM);
	printf("chi square test statistics=%lf\n",test_statistics);
	Free(C);
	Free(distances);
	Free(data->objects);
	Free(data);
	return 0;
}
