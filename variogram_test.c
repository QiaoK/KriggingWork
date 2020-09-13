#include "cluster.h"
#include "krigfunctions.h"
#include "random.h"
#include "clusterfunctions.h"
#include "datafunctions.h"

int main(){
	DTYPE angle_bound=M_PI/16, step_size=25,bound=step_size/2;
	SMOOTHING_TYPE smoothing_type=UNIFORM_SMOOTHING;
	DWORD steps=10,angle_steps=8;
	DWORD epochs=100,n_particles=200,variogram_type=ANISOTROPHY_POWER_VARIOGRAM;
	DWORD variogram_size=variogram_model_length(variogram_type);
	DTYPE c1=2,c2=2,alpha=0.729, gamma_hat;
	DWORD seed=555;
	DTYPE parameters[2];
	init_genrand(seed);
	DTYPE* C=Calloc(variogram_size,DTYPE);
	//Read SOCR data
	Objects* data=read_csv("SOCR_061708_NC_Data_Aquifer.csv",SPATIAL_DATA,FALSE);
	DWORD i;
	//Rescaling data
	for(i=0;i<data->size;i++){
		data->objects[i]->attribute*=1000;
	}
	printf("start\n");
	//Extract variogram samples.
	Samples* samples=variogram_sampling(data, bound, angle_bound, step_size, steps,angle_steps, C, smoothing_type);
	
	randomize_variogram(samples, C, variogram_type);
	C[3]=1.9;
	C[0]=14000;
	C=variogram_PSO(data,samples,epochs, n_particles, variogram_type, c1,c2,alpha);
	variogram_WLS(samples,C,5000, variogram_type,0.1);
	//C=variogram_PSO(data,samples,epochs, n_particles, variogram_type, c1,c2,alpha);
	//variogram_WLS_line_search(samples,C, 10, variogram_type);

	DWORD j;
	for(j=0;j<4;j++){
		printf("C[%lld]=%lf,",j,C[j]);
	}
	printf("\n");
	for(i=0;i<samples->size;i++){
		parameters[0]=samples->x[i];
		parameters[1]=samples->phi[i];
		gamma_hat=compute_variogram_by_parameters(parameters,C,variogram_type);
		//printf("x=%lf,y=%lf,f(x)=%lf,angle=%lf\n",samples->x[i],samples->y[i],gamma_hat,samples->phi[i]);
	}
	printf("real mse=%lf\n",evaluate_model(samples, C, variogram_type));
	C=Calloc(4,DTYPE);
	C[0]=14000;
	C[1]=38;
	C[2]=15;
	C[3]=1.99;
	for(i=0;i<samples->size;i++){
		parameters[0]=samples->x[i];
		parameters[1]=samples->phi[i];
		gamma_hat=compute_variogram_by_parameters(parameters,C,variogram_type);
		//printf("x=%lf,y=%lf,f(x)=%lf,angle=%lf\n",samples->x[i],samples->y[i],gamma_hat,samples->phi[i]);
	}
	printf("real mse=%lf\n",evaluate_model(samples, C, variogram_type));

	destroy_samples(samples);
	return 0;
}
