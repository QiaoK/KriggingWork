#include "cluster.h"
#include "krigfunctions.h"
#include "random.h"
#include "clusterfunctions.h"
#include "datafunctions.h"

int main(){
	//Read IGRA data
	Objects* data=read_IGRA("US_IGRA",60,221,TRUE);
	DWORD i,index=0;
	//Construct variogram
	VARIOGRAM_TYPE variogram_type=EXPONENTIAL_VARIOGRAM;
	DTYPE* C=Calloc(2,DTYPE);
	C[0]=54.84102;
	C[1]=153.71639;
	C[2]=0.05462576;
	DTYPE* distances=Calloc(2,DTYPE);
	distances[0]=30;
	distances[1]=DIS_UNCHECKED;
	Object** filtered=Calloc(data->size,Object*);
	for(i=0;i<data->size;i++){
		if(data->objects[i]->time==2014010100&&data->objects[i]->attribute!=-9999){
			filtered[index++]=data->objects[i];
			printf("x=%lf,y=%lf,t=%lf,time=%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->attribute,data->objects[i]->time);
		}
	}
	//printf("index=%lld\n",index);
	filtered=(Object**)realloc(filtered,sizeof(Object*)*index);
	Objects* objects=Calloc(1,Objects);
	objects->objects=filtered+1;
	objects->size=index-1;
	Matrix* weights1= krig_weights(objects,filtered[0],C,variogram_type);
	objects->size=10;
	Matrix* weights2= krig_weights(objects,filtered[0],C,variogram_type);
	DTYPE sum10=0;
	for(i=0;i<10;i++){
		sum10+=weights1->matrix[i][0];
	}
	printf("sum10=%lf\n",sum10);
	for(i=0;i<10;i++){
		printf("%lf,%lf,%lf\n",weights1->matrix[i][0]/sum10,weights2->matrix[i][0],distance(filtered[0]->spatial_coordinates,filtered[i+1]->spatial_coordinates));
	}
	return 0;
}
