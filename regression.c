#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "clustertype.h"
#include "cluster.h"
#include "krigfunctions.h"
#include "clusterfunctions.h"
#include "matrix.h"

Objects* getNearestNeighbors(Object* target,Objects* data,DWORD number){
	if(data->size<number||number<=0){
		return NULL;
	}
	DWORD i,j;
	DWORD max_index=0,empty=0;
	DTYPE max_dis=0,dis;
	Object** result=Calloc(number,Object*);
	for(i=0;i<number;i++){
		result[i]=data->objects[i];
		dis=distance(target->spatial_coordinates,data->objects[i]->spatial_coordinates);
		if(dis<=0){
			empty++;
		}
		if(dis>max_dis||dis<=0){
			max_dis=dis;
			max_index=i;
		}

	}
	//printf("max_dis=%lf,max_index=%lld,empty=%lld\n",max_dis,max_index,empty);
/*
	if(max_dis<=0){
		Free(result);
		return NULL;
	}
*/
	for(i=number;i<data->size;i++){
		dis=distance(target->spatial_coordinates,data->objects[i]->spatial_coordinates);
		if((max_dis<=0||dis<max_dis)&&dis>0){
			if(empty>0){
				empty--;
			}
			result[max_index]=data->objects[i];
			max_dis=dis;
			for(j=0;j<number;j++){
				if(max_dis<=0){
					break;
				}
				dis=distance(target->spatial_coordinates,result[j]->spatial_coordinates);
				if(dis>max_dis||dis<=0){
					max_dis=dis;
					max_index=j;
				}
			}
		}
	}
	if(empty>0){
		Free(result);
		return NULL;
	}
	Objects* wrapper=Calloc(1,Objects);
	wrapper->size=number;
	wrapper->objects=result;
	return wrapper;
}

Objects* locate_map(Objects* data,DTYPE time){
	Object** result=Calloc(data->size,Object*);
	DWORD i,index=0;
	for(i=0;i<data->size;i++){
		if(data->objects[i]->time==time){
			result[index++]=data->objects[i];
		}
	}
	result=(Object**)realloc(result,sizeof(Object*)*index);
	Objects* wrapper=Calloc(1,Objects);
	wrapper->size=index;
	wrapper->objects=result;
	return wrapper;
}

Objects* get_previous_time_stamps(Object* target,Objects* data,DWORD number){
	if(data->size<number||number<=0){
		return NULL;
	}
	DWORD i,j;
	DWORD max_index=0,empty=0;
	DTYPE max_dis=0,dis;
	Object** result=Calloc(number,Object*);
	for(i=0;i<number;i++){
		result[i]=data->objects[i];
		dis=target->time-data->objects[i]->time;
		if(dis<=0){
			empty++;
		}
		if(dis>max_dis||dis<=0){
			max_dis=dis;
			max_index=i;
		}

	}
/*
	if(max_dis<=0){
		Free(result);
		return NULL;
	}
*/
	for(i=number;i<data->size;i++){
		dis=target->time-data->objects[i]->time;
		if((max_dis<=0||dis<max_dis)&&dis>0){
			if(empty>0){
				empty--;
			}
			result[max_index]=data->objects[i];
			max_dis=dis;
			for(j=0;j<number;j++){
				if(max_dis<=0){
					break;
				}
				dis=target->time-result[j]->time;
				if(dis>max_dis||dis<=0){
					max_dis=dis;
					max_index=j;
				}
			}
		}
	}
	if(empty>0){
		Free(result);
		return NULL;
	}
	Objects* wrapper=Calloc(1,Objects);
	wrapper->size=number;
	wrapper->objects=result;
	return wrapper;
}
/*
 * Insertion sort. Assume the size of set being sorted is small.
 * Otherwise use qsort.
*/
void sort_by_distance(Object* target,Objects* objects){
	DWORD i,j,min_index;
	Object* temp;
	DTYPE min_dis;
	DTYPE dis;
	for(i=0;i<objects->size;i++){
		min_dis=INFINITY;
		min_index=i;
		for(j=i;j<objects->size;j++){
			dis=distance(target->spatial_coordinates,objects->objects[j]->spatial_coordinates);
			if(dis<min_dis){
				min_index=j;
				min_dis=dis;
			}
		}
		temp=objects->objects[i];
		objects->objects[i]=objects->objects[min_index];
		objects->objects[min_index]=temp;
	}
}

DTYPE predict_attribute(Objects* data,Object* target,DWORD time_lag,DWORD neighbor_number,Matrix* coefficients){
	DTYPE result=0;
	DWORD i,end=neighbor_number+time_lag;
	DTYPE dis;
	Objects *t_neighbors,*s_neighbors;
	s_neighbors=getNearestNeighbors(target,data,neighbor_number);
	if(s_neighbors==NULL&&neighbor_number>0){
		return INFINITY;
	}
	//printf("check 1\n");
	t_neighbors=get_previous_time_stamps(target,data,time_lag);
	if(t_neighbors==NULL&&time_lag>0){
		if(s_neighbors!=NULL){
			Free(s_neighbors->objects);
			Free(s_neighbors);
		}
		return INFINITY;
	}
	if(s_neighbors!=NULL){
		sort_by_distance(target,s_neighbors);
	}
	for(i=0;i<time_lag;i++){
		result+=coefficients->matrix[i][0]*t_neighbors->objects[i]->attribute;
	}
	for(i=time_lag;i<end;i++){
		dis=distance(s_neighbors->objects[i]->spatial_coordinates,target->spatial_coordinates);
		result+=coefficients->matrix[i][0]*s_neighbors->objects[i]->attribute/dis;
	}
	return result+coefficients->matrix[end][0];
}

Matrix* construct_st_model(Objects* data,DWORD time_lag,DWORD neighbor_number){
	DWORD i,j,end=neighbor_number+time_lag;
	Matrix* f=create_matrix(end+1,1);
	Matrix* ff=create_matrix(end+1,end+1);
	Matrix* fy=create_matrix(end+1,1);
	Matrix* temp;
	DTYPE dis;
	Objects *t_neighbors,*s_neighbors;
	for(i=0;i<ff->n_row;i++){
		fy->matrix[i][0]=0;
		for(j=0;j<ff->n_column;j++){
			ff->matrix[i][j]=0;
		}
	}
	f->matrix[end][0]=1;
	for(i=0;i<data->size;i++){
		s_neighbors=getNearestNeighbors(data->objects[i],data,neighbor_number);
		if(s_neighbors==NULL&&neighbor_number>0){
			continue;
		}
		//printf("check 1\n");
		t_neighbors=get_previous_time_stamps(data->objects[i],data,time_lag);
		if(t_neighbors==NULL&&time_lag>0){
			if(s_neighbors!=NULL){
				Free(s_neighbors->objects);
				Free(s_neighbors);
			}
			continue;
		}
		if(s_neighbors!=NULL){
			sort_by_distance(data->objects[i],s_neighbors);
		}
		for(j=0;j<time_lag;j++){
			f->matrix[j][0]=t_neighbors->objects[j]->attribute;
			fy->matrix[j][0]+=data->objects[i]->attribute*f->matrix[j][0];
		}
		for(j=time_lag;j<end;j++){
			dis=distance(s_neighbors->objects[j]->spatial_coordinates,data->objects[i]->spatial_coordinates);
			f->matrix[j][0]=s_neighbors->objects[j]->attribute/dis;
			fy->matrix[j][0]+=data->objects[i]->attribute*f->matrix[j][0];
		}
		fy->matrix[end][0]+=data->objects[i]->attribute;
		temp=vector_self_multiplication(f);
		add_to_matrix(ff,temp);
		//printf("check 3\n");
		destroy_matrix(temp);
		if(t_neighbors!=NULL){
			Free(t_neighbors->objects);
			Free(t_neighbors);
		}
		if(s_neighbors!=NULL){
			Free(s_neighbors->objects);
			Free(s_neighbors);
		}
	}
	printf("ff\n");
	print_matrix(ff);
	printf("fy\n");
	print_matrix(fy);
	//printf("check\n");
	Matrix* coefficients=solve_linear_system(ff,fy);
	//print_matrix(matrix_multiplication(ff,coefficients));
	return coefficients;
}
/*
void main(){
	Cluster* cluster=create_cluster();
	DWORD i;
	for(i=0;i<10;i++){
		Object* object=Calloc(1,Object);
		object->time=i;
		object->spatial_coordinates=Calloc(2,DTYPE);
		object->spatial_coordinates[0]=i;
		object->spatial_coordinates[1]=i;
		add_to_cluster(cluster,object);
	}
	Objects* data=cluster_to_objects(cluster);
	Objects* map=locate_map(data,2);
	printf("size of map=%ld\n",map->size);
	Objects* neighbors=getNearestNeighbors(data->objects[6],data,4);
	sort_by_distance(data->objects[6],neighbors);
	if(neighbors==NULL){
		printf("insufficient data\n");
	}else{
		for(i=0;i<neighbors->size;i++){
			printf("time=%lf\n",neighbors->objects[i]->time);
		}
	}
	Free(map->objects);
	Free(map);
	destroy_cluster(cluster);
}
*/
