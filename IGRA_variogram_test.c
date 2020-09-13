/*
*  Copyright (C) 2017, Northwestern University.
*/
#include "cluster.h"
#include "krigfunctions.h"
#include "random.h"
#include "clusterfunctions.h"
#include "datafunctions.h"

int main(void){
	DTYPE angle_bound=M_PI*2+1, step_size=3,bound=step_size/2;
	SMOOTHING_TYPE smoothing_type=UNIFORM_SMOOTHING;
	DWORD steps=10,angle_steps=1;
	DWORD epochs=100,n_particles=200,variogram_type=SPHERICAL_VARIOGRAM;
	Objects* data=read_IGRA("US_IGRA",60,221,TRUE);
	DWORD i,j,index=0;
	DTYPE* C=Calloc(3,DTYPE);
	DTYPE increment=-88;
	DWORD variogram_size=variogram_model_length(variogram_type);
	DTYPE c1=2,c2=2,alpha=0.729,gamma_hat;
	DWORD seed=5555;
	init_genrand(seed);

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
		Objects* filtered_data=Calloc(1,Objects);
		filtered_data->size=index;
		filtered_data->objects=filtered;
		Samples* samples=variogram_sampling(filtered_data, bound, angle_bound, step_size, steps,angle_steps, C, smoothing_type);
		//print_samples(samples);
		randomize_variogram(samples, C, variogram_type);
		
		for(i=0;i<3;i++){
			printf("C[%lld]=%lf,",i,C[i]);
		}
		printf("real mse=%lf\n",evaluate_model(samples, C, variogram_type));
		
		//C=variogram_PSO(data,samples,epochs, n_particles, variogram_type, c1,c2,alpha);
		variogram_WLS(samples,C,2000, variogram_type,.5);

		for(i=0;i<3;i++){
			printf("C[%lld]=%lf,",i,C[i]);
		}
		printf("real mse=%lf\n",evaluate_model(samples, C, variogram_type));
		//write_variogram_result("test.txt", TRUE, samples, C, variogram_type);
		destroy_samples(samples);
		Free(filtered);
		Free(filtered_data);
	}
	Free(C);
}
