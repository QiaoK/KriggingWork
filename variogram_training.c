#include "matrix.h"
#include "clusterfunctions.h"
#include "krigfunctions.h"
#include "random.h"
#define VARIOGRAM_EXP_BOUND .05

Samples* create_samples(DWORD size){
	//Malloc space for sample
	Samples* samples=Calloc(1,Samples);
	//Initialize its components
	samples->x=Calloc(size,DTYPE);
	samples->y=Calloc(size,DTYPE);
	samples->phi=Calloc(size,DTYPE);
	samples->N=Calloc(size,DWORD);
	samples->size=size;
	return samples;
}

void destroy_samples(Samples *samples){
	Free(samples->x);
	Free(samples->y);
	Free(samples->phi);
	Free(samples->N);
	Free(samples);
}

void print_samples(Samples* samples){
	DWORD i;
	for(i=0;i<samples->size;i++){
		printf("x=%lf,y=%lf,angle=%lf,N=%lld\n",samples->x[i],samples->y[i],samples->phi[i],samples->N[i]);
	}
}

DTYPE normal_kernel_smoothing(DTYPE value,DTYPE *C){
	return exp(-value*value/(2*C[0]*C[0]));
}

DTYPE kernel_smoother(DTYPE value,DTYPE *C,SMOOTHING_TYPE smoothing_type){
	switch(smoothing_type){
		case GAUSSIAN_SMOOTHING:
			return normal_kernel_smoothing(value,C);
		case UNIFORM_SMOOTHING:
			return 1;
		default:
			return 1;
	}
}

void smooth_variance_estimation(Samples* samples,DWORD index,Objects* objects,DTYPE bound,DTYPE angle_bound,DTYPE lag,DTYPE angle,DTYPE *C,SMOOTHING_TYPE smoothing_type){
	DWORD i,j;
	DTYPE result=0,dis,arc_dis,temp1,temp2,smoother=0;
	samples->N[index]=0;
	for(i=0;i<objects->size;i++){
		for(j=i+1;j<objects->size;j++){
			dis=distance(objects->objects[i]->spatial_coordinates,objects->objects[j]->spatial_coordinates);
			arc_dis=atan((objects->objects[i]->spatial_coordinates[1]-objects->objects[j]->spatial_coordinates[1])/(objects->objects[i]->spatial_coordinates[0]-objects->objects[j]->spatial_coordinates[0]));
			temp1=objects->objects[i]->attribute-objects->objects[j]->attribute;
			if(fabs(dis-lag)<bound&&fabs(arc_dis-angle)<angle_bound){
				//printf("arc_dis=%lf,angle=%lf,angle_bound=%lf,dis=%lf,lag=%lf\n",arc_dis,angle,angle_bound,dis,lag);
				temp2=kernel_smoother(temp1,C,smoothing_type);
				result+=temp1*temp1*temp2;
				smoother+=temp2;
				samples->N[index]++;
			}

		}
	}
	if(smoother==0){
		samples->y[index]=-1;
	}else{
		samples->y[index]=result/smoother;
	}
	//printf("smoother=%lf\n",smoother);
}

Samples* variogram_sampling(Objects *objects,DTYPE bound,DTYPE angle_bound,DTYPE step_size,DWORD steps,DWORD angle_steps,DTYPE *C,SMOOTHING_TYPE smoothing_type){
	Samples* samples=create_samples(steps*angle_steps);
	DWORD i,j,index=0;
	DTYPE lag=step_size/2,angle;
	for(i=0;i<steps;i++){
		angle=0;
		for(j=0;j<angle_steps;j++){
			samples->x[index]=lag;
			samples->phi[index]=angle;
			smooth_variance_estimation(samples,index,objects,bound,angle_bound,lag,angle,C,smoothing_type);
			angle+=M_PI/angle_steps;
			if(samples->y[index]>0){
				index++;
			}
		}
		lag+=step_size;
	}
	samples->x=(DTYPE*)realloc(samples->x,sizeof(DTYPE)*index);
	samples->y=(DTYPE*)realloc(samples->y,sizeof(DTYPE)*index);
	samples->phi=(DTYPE*)realloc(samples->phi,sizeof(DTYPE)*index);
	samples->size=index;
	return samples;
}

void randomize_variogram(Samples* samples, DTYPE* C, VARIOGRAM_TYPE variogram_type){
	DTYPE min=samples->y[0],max=samples->y[0],mean=samples->y[0];
	DWORD i;
	for(i=1;i<samples->size;i++){
		if(samples->y[i]>max){
			max=samples->y[i];
		}
		if(samples->y[i]<min){
			min=samples->y[i];
		}
		mean+=samples->y[i];
	}
	mean/=samples->size;
	
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			C[0]=(genrand_real2()*.25+.75)*min;
			C[1]=genrand_real2()*max+mean*.5;
			C[2]=-(1+genrand_real2())*log(VARIOGRAM_EXP_BOUND)/mean;
			break;
		}
		case SPHERICAL_VARIOGRAM:{
			C[0]=(genrand_real2()*.25+.75)*min;
			C[1]=sqrt(genrand_real2()*max);
			C[2]=(1+genrand_real2())*samples->x[samples->size-1];
			break;
		}
		case POWER_VARIOGRAM:{
			
			break;
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{
			C[0]=min+genrand_real2()*(mean-min)/2;
			C[1]=genrand_real2()*100;
			C[2]=genrand_real2()*100;
			C[3]=1.5+genrand_real2()/2;
			break;
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			
			break;
		}

		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
			
			break;
		}
		default:{
			
			break;
		}
	}
}

DWORD variogram_model_length(VARIOGRAM_TYPE variogram_type){
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			return 3;
			break;
		}
		case SPHERICAL_VARIOGRAM:{
			return 3;
			break;
		}
		case POWER_VARIOGRAM:{
			break;
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{
			return 4;
			break;
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			
			break;
		}
		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
	
			break;
		}
		default:{
			
			break;
		}
	}
	return -1;
}

Particle* create_variogram_particle(DWORD size){
	Particle* particle=Calloc(1,Particle);
	particle->C=Calloc(size,DTYPE);
	particle->local_best=Calloc(size,DTYPE);
	memset(particle->local_best,0,size*sizeof(DTYPE));
	particle->velocity=Calloc(size,DTYPE);
	memset(particle->velocity,0,size*sizeof(DTYPE));
	return particle;
}

Particle** create_variogram_particles(DWORD size,VARIOGRAM_TYPE variogram_size){
	Particle **particles=Calloc(size,Particle*);
	DWORD i;
	for(i=0;i<size;i++){
		particles[i]=create_variogram_particle(variogram_size);
	}
	return particles;
}

void destroy_variogram_particle(Particle* particle){
	Free(particle->C);
	Free(particle->local_best);
	Free(particle->velocity);
	Free(particle);
}

void destroy_variogram_particles(Particle** particles,DWORD size){
	DWORD i;
	for(i=0;i<size;i++){
		destroy_variogram_particle(particles[i]);
	}
	Free(particles);
}

DTYPE evaluate_model(Samples* samples, DTYPE* C, VARIOGRAM_TYPE variogram_type){
	DWORD i;
	DTYPE parameters[2];
	DTYPE gamma_hat,result=0,temp;
	for(i=0;i<samples->size;i++){
		parameters[0]=samples->x[i];
		parameters[1]=samples->phi[i];
		gamma_hat=compute_variogram_by_parameters(parameters,C,variogram_type);
		temp=gamma_hat/samples->y[i]-1;
		result+=samples->N[i]*temp*temp;
		//printf("gamma_hat=%lf\n",samples->N[i]*temp*temp);
	}
	return result;
}

void compute_variogram_gradient(Samples* samples,VARIOGRAM_TYPE variogram_type,DTYPE*C,DTYPE* gradients){
	DWORD i;
	DTYPE estimation;
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			DTYPE base;
			DWORD count=0;
			gradients[0]=0;
			gradients[1]=0;
			gradients[2]=0;
			DTYPE parameters[1];
			for(i=0;i<samples->size;i++){
				parameters[0]=samples->x[i];
				estimation=compute_variogram_by_parameters(parameters,C,variogram_type);
				base=2*samples->N[i]*(1-estimation/samples->y[i])/samples->y[i];
				//gradients[1]+=base*2*C[1]*(1-exp(-C[2]*samples->x[i]));
				gradients[2]+=base*C[1]*C[1]*C[2]*samples->x[i]*exp(-C[2]*samples->x[i]);
				count+=samples->N[i];
			}
			gradients[0]/=count;
			gradients[1]/=count;
			gradients[2]/=count;
			break;
		}
		case SPHERICAL_VARIOGRAM:{
			DTYPE base;
			DWORD count=0;
			gradients[0]=0;
			gradients[1]=0;
			gradients[2]=0;
			DTYPE parameters[1];
			for(i=0;i<samples->size;i++){
				parameters[0]=samples->x[i];
				estimation=compute_variogram_by_parameters(parameters,C,variogram_type);
				base=2*samples->N[i]*(1-estimation/samples->y[i])/samples->y[i];
				gradients[1]+=base*2*C[1]*(1.5*samples->x[i]/C[2]-0.5*samples->x[i]*samples->x[i]*samples->x[i]/(C[2]*C[2]*C[2]));
				gradients[2]+=base*C[1]*C[1]*1.5*(samples->x[i]*samples->x[i]*samples->x[i]/(C[2]*C[2]*C[2]*C[2])-samples->x[i]/(C[2]*C[2]));
				count+=samples->N[i];
			}
			gradients[0]/=count;
			gradients[1]/=count;
			gradients[2]/=count;
			break;
		}
		case POWER_VARIOGRAM:{
			
			break;
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{

			DTYPE parameters[2];
			DTYPE base;
			DTYPE r2,phi;
			DWORD count=0;
			gradients[0]=0;
			gradients[1]=0;
			gradients[2]=0;
			gradients[3]=0;
			for(i=0;i<samples->size;i++){
				parameters[0]=samples->x[i];
				parameters[1]=samples->phi[i];
				r2=samples->x[i]*samples->x[i];
				phi=samples->phi[i];
				estimation=compute_variogram_by_parameters(parameters,C,variogram_type);
				base=-2*samples->N[i]*(1-estimation/samples->y[i])/samples->y[i];
				//gradients[0]+=base;
				gradients[1]+=base*2*r2*cos(M_PI/4-phi)*pow(C[1],2/C[3]-1)/C[3];
				gradients[2]+=base*2*r2*cos(M_PI/4+phi)*pow(C[2],2/C[3]-1)/C[3];
				//gradients[3]+=base*(-r2*(pow(C[1],2/C[3])*log(C[1])*cos(M_PI/4-phi)+pow(C[2],2/C[3])*log(C[2])*cos(M_PI/4+phi)))/(C[3]*C[3]);
				count+=samples->N[i];
			}
			
			gradients[0]/=count;
			gradients[1]/=count;
			gradients[2]/=count;
			gradients[3]/=count;
			
			//printf("gradient[0]=%lf,gradient[1]=%lf,gradient[2]=%lf,gradient[3]=%lf\n",gradients[0],gradients[1],gradients[2],gradients[3]);
			break;
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			
			break;
		}

		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
			
			break;
		}
		default:{
			
			break;
		}
	}
}

void variogram_WLS(Samples* samples,DTYPE* C,DWORD epochs, VARIOGRAM_TYPE variogram_type,DTYPE learning_rate){
	DWORD i,j;
	DWORD variogram_size=variogram_model_length(variogram_type);
	DTYPE* gradients=Calloc(variogram_size,DTYPE);

	for(i=0;i<epochs;i++){
		compute_variogram_gradient(samples,variogram_type,C,gradients);
		for(j=0;j<variogram_size;j++){
			C[j]+=gradients[j]*learning_rate;
		}
		/*
		for(j=0;j<3;j++){
			printf("gradients[%lld]=%lf,",j,gradients[j]);
		}
		printf("\n");
		*/
	}
	Free(gradients);
}

DTYPE line_search_alpha(Samples* samples,DTYPE* C, DTYPE* d, DTYPE* gradients, VARIOGRAM_TYPE variogram_type){
	DWORD i;
	DTYPE alpha=0;
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			break;
		}
		case POWER_VARIOGRAM:{
			
			break;
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{
			DTYPE s_derivative1=0;
			DTYPE s_derivative2=0;
			DTYPE estimation;
			DTYPE parameters[2];
			DTYPE base;
			DTYPE r2,phi;
			DWORD count=0;
			for(i=0;i<samples->size;i++){
				parameters[0]=samples->x[i];
				parameters[1]=samples->phi[i];
				r2=samples->x[i]*samples->x[i];
				phi=samples->phi[i];
				estimation=compute_variogram_by_parameters(parameters,C,variogram_type);
				base=-2*samples->N[i]*(1-samples->y[i]/estimation)*samples->y[i]/(estimation*estimation);
				//printf("base=%lf,%lf\n",base,(C[3]/2-1)/C[3]);
				s_derivative1+=base*2*r2*cos(M_PI/4-phi)*pow(C[1],2/C[3]-2)*(2/C[3]-1)/C[3];
				s_derivative2+=base*2*r2*cos(M_PI/4+phi)*pow(C[2],2/C[3]-2)*(2/C[3]-1)/C[3];
				count+=samples->N[i];
			}
			s_derivative1/=count;
			s_derivative2/=count;
			DTYPE denominator=d[1]*d[1]*s_derivative1+d[2]*d[2]*s_derivative2;
			DTYPE nominator=d[1]*(s_derivative1*C[1]+gradients[1])+d[2]*(s_derivative2*C[2]+gradients[2]);
			printf("alpha,d[1]=%lf,d[2]=%lf,s_derivative1=%lf,s_derivative2=%lf,gradients[1]=%lf,gradients[2]=%lf,denominator=%lf,nominator=%lf\n",d[1],d[2],s_derivative1,s_derivative2,gradients[1],gradients[2],denominator,nominator);
			alpha=-nominator/denominator;
			break;
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			
			break;
		}

		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
			
			break;
		}
		default:{
			
			break;
		}
	}
	return alpha;
}

DTYPE line_search_beta(Samples* samples,DTYPE* C, DTYPE* d, DTYPE* gradients, VARIOGRAM_TYPE variogram_type){
	DWORD i;
	DTYPE beta=0;
	switch(variogram_type){
		case EXPONENTIAL_VARIOGRAM:{
			break;
		}
		case POWER_VARIOGRAM:{
			
			break;
		}
		case ANISOTROPHY_POWER_VARIOGRAM:{
			DTYPE s_derivative1=0;
			DTYPE s_derivative2=0;
			DTYPE estimation;
			DTYPE parameters[2];
			DTYPE base;
			DTYPE r2,phi;
			DWORD count=0;
			for(i=0;i<samples->size;i++){
				parameters[0]=samples->x[i];
				parameters[1]=samples->phi[i];
				r2=samples->x[i]*samples->x[i];
				phi=samples->phi[i];
				estimation=compute_variogram_by_parameters(parameters,C,variogram_type);
				base=-2*samples->N[i]*(1-samples->y[i]/estimation)*samples->y[i]/(estimation*estimation);
				//printf("base=%lf,%lf\n",base,(C[3]/2-1)/C[3]);
				s_derivative1+=base*2*r2*cos(M_PI/4-phi)*pow(C[1],2/C[3]-2)*(2/C[3]-1)/C[3];
				s_derivative2+=base*2*r2*cos(M_PI/4+phi)*pow(C[2],2/C[3]-2)*(2/C[3]-1)/C[3];
				count+=samples->N[i];
			}
			s_derivative1/=count;
			s_derivative2/=count;
			DTYPE denominator=d[1]*d[1]*s_derivative1+d[2]*d[2]*s_derivative2;
			DTYPE nominator=d[1]*gradients[1]*s_derivative1+d[2]*gradients[2]*s_derivative2;
			printf("beta,d[1]=%lf,d[2]=%lf,s_derivative1=%lf,s_derivative2=%lf,gradients[1]=%lf,gradients[2]=%lf,denominator=%lf,nominator=%lf\n",d[1],d[2],s_derivative1,s_derivative2,gradients[1],gradients[2],denominator,nominator);
			beta=-nominator/denominator;
			break;
		}
		case ST_SPHERICAL_PRODUCT_VARIOGRAM:{
			
			break;
		}

		case ST_EXPONENTIAL_PRODUCT_VARIOGRAM:{
			
			break;
		}
		default:{
			
			break;
		}
	}
	return beta;
}

void variogram_WLS_line_search(Samples* samples,DTYPE *C, DWORD epochs, VARIOGRAM_TYPE variogram_type){
	DWORD i,j;
	DWORD variogram_size=variogram_model_length(variogram_type);
	DTYPE* gradients=Calloc(variogram_size,DTYPE);
	DTYPE* d=Calloc(variogram_size,DTYPE);
	DTYPE alpha,beta;
	for(i=0;i<epochs;i++){
		compute_variogram_gradient(samples,variogram_type,C,gradients);
		if(i==0){
			for(j=0;j<variogram_size;j++){
				d[j]=-gradients[j];
			}
		}else{
			beta=line_search_beta(samples,C, d, gradients, variogram_type);
			printf("beta=%lf\n",beta);
			for(j=0;j<variogram_size;j++){
				d[j]=-gradients[j]+beta*d[j];
			}
		}
		alpha=line_search_alpha(samples,C, d, gradients, variogram_type);
		printf("alpha=%lf\n",alpha);
		for(j=0;j<variogram_size;j++){
			C[j]+=alpha*d[j];
		}
		for(j=0;j<4;j++){
			printf("C[%lld]=%lf,",j,C[j]);
		}
		printf("\n");
	}
	Free(gradients);
	Free(d);
}

void copy_particle(DTYPE* src,DTYPE* dest,DTYPE variogram_size){
	memcpy(dest,src,variogram_size*sizeof(DTYPE));
}

void update_particle(Particle* particle,DTYPE* global_best,DTYPE c1,DTYPE c2,DTYPE alpha,DTYPE phi1,DTYPE phi2,DWORD variogram_size){
	DWORD i;
	DTYPE l_best,g_best,loc;
	for(i=1;i<variogram_size;i++){
		l_best=particle->local_best[i];
		g_best=global_best[i];
		loc=particle->C[i];
		particle->velocity[i]=loc*alpha+phi1*c1*(l_best-loc)+phi2*c2*(g_best-loc);
		particle->C[i]+=particle->velocity[i];
	}
}

DTYPE* variogram_PSO(Objects *objects,Samples* samples,DWORD epochs, DWORD n_particles, VARIOGRAM_TYPE variogram_type, DTYPE c1,DTYPE c2,DTYPE alpha){
	DWORD i,j;
	DTYPE *temp_best=NULL,*global_best;
	DTYPE best=-1,temp;
	DTYPE *phi1=Calloc(n_particles,DTYPE);
	DTYPE *phi2=Calloc(n_particles,DTYPE);
	DWORD variogram_size=variogram_model_length(variogram_type);
	DTYPE* distances=Calloc(2,DTYPE);
	distances[0]=DIS_UNCHECKED;
	distances[1]=DIS_UNCHECKED;
	//Create cluster and insert points
        Cluster *cluster=create_cluster();
	for(i=0;i<objects->size;i++){
		add_to_cluster(cluster,objects->objects[i]);
	}
	//initialize particles
	Particle **particles=create_variogram_particles(n_particles, variogram_size);
	for(i=0;i<n_particles;i++){
		randomize_variogram(samples, particles[i]->C, variogram_type);
		particles[i]->C[0]=25;
	}
	//compute g_best
	for(i=0;i<n_particles;i++){
		particles[i]->mse=evaluate_model(samples, particles[i]->C, variogram_type);
		for(j=0;j<variogram_size;j++){
			//printf("C[%lld]=%lf,",j,particles[i]->C[j]);
		}
		//printf("mse=%lf\n",particles[i]->mse);
	}
	for(i=0;i<n_particles;i++){
		if(best==-1||particles[i]->mse<best){
			best=particles[i]->mse;
			temp_best=particles[i]->C;
		}
	}
	global_best=Calloc(variogram_size,DTYPE);
	copy_particle(temp_best,global_best,variogram_size);
	//PSO algorithm
	for(i=0;i<epochs;i++){
                //DTYPE cluster_error=sqrt(sum_krig_variance(cluster,global_best,distances,variogram_type));
		//printf("epoch %lld, g_best=%lf,Kriging error=%lf\n",i,best,cluster_error);
		for(j=0;j<n_particles;j++){
			phi1[j]=genrand_real2();
			phi2[j]=genrand_real2();
		}
		temp_best=NULL;
		for(j=0;j<n_particles;j++){
			update_particle(particles[j],global_best,c1,c2,alpha,phi1[j],phi2[j],variogram_size);
			temp=evaluate_model(samples, particles[j]->C, variogram_type);
			if(temp<particles[j]->mse){
				particles[j]->mse=temp;
				copy_particle(particles[j]->C,particles[j]->local_best,variogram_size);
				if(particles[j]->mse<best){
					best=particles[j]->mse;
					temp_best=particles[j]->C;
				}
			}
		}
		if(temp_best!=NULL){
			copy_particle(temp_best,global_best,variogram_size);
		}

	}
	destroy_variogram_particles(particles,n_particles);
	Free(phi1);
	Free(phi2);
	Free(distances);
	return global_best;
}

