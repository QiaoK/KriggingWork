/*
 * Copyright (C) 2016, Northwestern University.
 * This file contains the core Kriging clustering algorithm.
 * See also krigfunctions.h
*/

#include "clusterfunctions.h"
#include "krigfunctions.h"

Clusters* krig_clustering(Object** data,DWORD size,DTYPE bound,DTYPE *C,DTYPE* max_distance,VARIOGRAM_TYPE variogram_type){
	/*Initialize variables*/
	struct timeval start,end; /*Timing variables*/
	DWORD i,j,k=1,counter,change; /*Temporary variables*/
	DWORD revision_time=0,filter_time=0; /*Counter variables for iterations.*/
	//Cluster* clone;
	Node *temp,*current,*previous; /*Temporary variables for iterating linked list (cluster)*/
	DTYPE predict;
	Cluster** clusters=Calloc(size,Cluster*);
	clusters[0]=create_cluster();
	DWORD filter_steps=0;
	DWORD revision_steps=0;
	add_to_cluster(clusters[0],data[0]);
	/*Add all data to a single cluster.*/
	for(i=1;i<size;i++){
		add_to_cluster(clusters[0],data[i]);
		clusters[i]=NULL;
	}
	for(i=0;i<k;i++){
		change=TRUE;
		//printf("Begin to filter--------------------------------\n");
		/*Filtering phase starts.*/
		/*For each clusters that are not filtered yet.*/
		while(change&&clusters[i]->size>1){
			current=clusters[i]->head;
			previous=NULL;
			change=FALSE;
			Node* tail=clusters[i]->tail;
			//clone=clone_cluster(clusters[i]);
			counter=0;
			/*For each element in the cluster*/
			while(clusters[i]->size>1&&current!=NULL){
				filter_steps++;
				counter++;
				gettimeofday(&start,NULL);
				//printf("cluster=%lld,size=%lld\n",i,clusters[i]->size);
				/*Remove the element from the cluster*/
				remove_from_cluster(clusters[i],previous,current);
				/*Use the rest of points in the same cluster to predict its non-spatial attribute*/
				predict=krig_normalize(clusters[i],current->object,C,max_distance,variogram_type);
				//predict=krig_normalize(clone,current->object,C,max_distance,variogram_type);
				//printf("cluster size=%lld,filter step=%lld,normalized error=%lf\n",clusters[i]->size,counter,predict);
				//printf("var=%lf,predicted=%lf,real=%lf\n",krig_variance(clusters[i],current->object,C,max_distance,variogram_type),krig_prediction(clusters[i],current->object,C,max_distance,variogram_type),current->object->attribute);
				/*Insert the point being filtered back to the cluster.*/
				insert_to_cluster(clusters[i],previous,current);
				if(clusters[i]->tail!=tail){
					/*Should not reach here, testing purpose*/
					printf("warning\n");
				}
				if(fabs(predict)>bound){
					/*Filtered out case*/
					change=TRUE;
					/*Create a new cluster if the algorithm has not*/
					if(clusters[i+1]==NULL){
						k++;
						clusters[i+1]=create_cluster();
					}
					/*Add the filtered point to next cluster*/
					add_to_cluster(clusters[i+1],current->object);
					current->object->neighbors=-1;
					remove_from_cluster(clusters[i],previous,current);
					temp=current;
					current=current->next;
					Free(temp);
				}else{
					/*The point stays in the cluster by passing the filtering phase*/
					previous=current;
					current=current->next;
				}
				/*Timing for filtering phase*/
				gettimeofday(&end,NULL);
				filter_time+=(end.tv_sec-start.tv_sec);

			}
			//destroy_cluster(clone);
			//filter_cluster(clusters[i],clusters[i+1]);
		}
		//printf("i=%lld,size=%lld\n",i,clusters[i]->size);'
		/*If nothing has been filtered, exit the loop*/
		if(clusters[i+1]==NULL||clusters[i+1]->size==1){
			break;
		}
		//printf("Begin to revise--------------------------------\n");
		/*Reinforcement (revising) phase starts*/
		change=TRUE;
		while(change&&clusters[i+1]->size>1){
			change=FALSE;
			current=clusters[i+1]->head;
			previous=NULL;
			counter=0;
			//Iterate all points in the cluster that contains points filtered out at filtering phase.
			while(current!=NULL){
				counter++;
				revision_steps++;
				gettimeofday(&start,NULL);
				/*Use all points in the cluster that was filtered to predict its non-spatial attribute*/
				predict=krig_normalize(clusters[i],current->object,C,max_distance,variogram_type);
				//printf("cluster size=%lld,predict for revising=%lf,step=%lld\n",clusters[i+1]->size,predict,counter);
				//remove_from_cluster(clusters[i+1],previous,current);
				add_to_cluster_front(clusters[i],current->object);
				//printf("check5\n");
				if(fabs(predict)>bound||!krig_consistency(clusters[i],bound,C,max_distance,variogram_type)){
					/*If not consistent, keep the point in the current cluster*/
					remove_cluster_front(clusters[i]);
					previous=current;
					Objects* objects=get_adjacent_objects(clusters[i],current->object,max_distance);
					for(j=0;j<objects->size;j++){
						objects->objects[j]->neighbors=-1;
					}
					Free(objects->objects);
					Free(objects);
					current->object->neighbors=-1;
					current=current->next;
				}else{
					/*If consistent, put the point back to the cluster that it was filtered out from.*/
					change=TRUE;
					remove_from_cluster(clusters[i+1],previous,current);
					current->object->neighbors=-1;
					current=current->next;
				}
				gettimeofday(&end,NULL);
				revision_time+=(end.tv_sec-start.tv_sec);
			}
			//printf("%lld\n",clusters[i+1]->size);
		}
		/*If there is only one point left in the reinforced cluster after reinforcement phase, algorithm exits*/
		if(clusters[i+1]->size<=1){
			break;
		}
		//printf("i=%lld,size=%lld\n",i,clusters[i]->size);
	}
	/*Point out timing statistics*/
	printf("filters=%lld,revisions=%lld,filter time=%lld s,revision time =%lld s\n",filter_steps,revision_steps,filter_time,revision_time);
	clusters=(Cluster**)realloc(clusters,sizeof(Cluster*)*k);
	Clusters *result=Calloc(1,Clusters);
	result->size=k;
	result->clusters=clusters;
	return result;
}
