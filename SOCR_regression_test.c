#include "matrix.h"
#include "krigfunctions.h"
#include "random.h"
#include "clusterfunctions.h"
#include "datafunctions.h"
#include "regression.h"

void print_data(Objects* data){
	DWORD i;
	for(i=0;i<data->size;i++){
		printf("x=%lf,y=%lf,attribute=%lf\n",data->objects[i]->spatial_coordinates[0],data->objects[i]->spatial_coordinates[1],data->objects[i]->attribute);
	}
}

int main(void){
	Objects* data=read_csv("SOCR_061708_NC_Data_Aquifer.csv",SPATIAL_DATA,FALSE);
	DWORD i;
	for(i=0;i<data->size;i++){
		data->objects[i]->attribute*=1000;
	}
	//Objects* neighbors=getNearestNeighbors(data->objects[6],data,4);
	shuffle_objects(data->objects,data->size);
	DWORD size=data->size;
	data->size=data->size*9/10;
	Matrix* coefficients=construct_st_model(data,0,7);
	data->size=size;
	print_matrix(coefficients);
	for(i=data->size*9/10;i<data->size;i++){
		printf("%lf,%lf\n",data->objects[i]->attribute,predict_attribute(data,data->objects[i],0,7,coefficients));
	}
	destroy_matrix(coefficients);
	return 0;
}
