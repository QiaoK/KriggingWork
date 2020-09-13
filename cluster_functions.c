/*
*  Copyright (C) 2016, Northwestern University.
*/
#include "clusterfunctions.h"


Cluster* create_cluster(){
	//Malloc space for this cluster
	Cluster* cluster=Calloc(1,Cluster);
	//Initialize its components to represent empty cluster
	cluster->size=0;
	cluster->tail=NULL;
	cluster->head=NULL;
	return cluster;
}

void add_to_cluster(Cluster* cluster,Object* object){
	if(cluster->size==0){
		//If the cluster is empty, insert the element as the first element.
		cluster->tail=Calloc(1,Node);
		cluster->tail->object=object;
		cluster->size=1;
		cluster->head=cluster->tail;
		cluster->tail->next=NULL;
		cluster->sum=object->attribute;
	}else{
		//If the cluster is not empty, insert the element to the back of this cluster (link list).
		Node* tail=Calloc(1,Node);
		tail->next=NULL;
		cluster->tail->next=tail;
		cluster->tail=tail;
		cluster->tail->object=object;
		cluster->size+=1;
		cluster->sum+=object->attribute;
	}
}

Cluster* clone_cluster(Cluster* cluster){
	//Initialize a new cluster.
	Cluster* result=create_cluster();
	//Copy everything from the input to the new cluster.
	Node* front=cluster->head;
	while(front!=NULL){
		add_to_cluster(result,front->object);
		front=front->next;
	}
	return result;
}

Object** get_objects(Cluster* cluster){
	//create an array of objects
	Object** objects=Calloc(cluster->size,Object*);
	//fill the array with elements in the cluster.
	Node* front=cluster->head;
	DWORD i;
	for(i=0;i<cluster->size;i++){
		objects[i]=front->object;
		front=front->next;
	}
	return objects;
}

Objects* cluster_to_objects(Cluster* cluster){
	//encapsulate an array of objects with size.
	Objects* objects=Calloc(1,Objects);
	objects->objects=get_objects(cluster);
	objects->size=cluster->size;
	return objects;
}

void destroy_cluster(Cluster* cluster){
	DWORD i;
	Node *next=cluster->head;
	Node *temp;
	//Free individual node, which is the container for object pointer, in the cluster 
	for(i=0;i<cluster->size;i++){
		temp=next->next;
		Free(next);
		next=temp;
	}
	Free(cluster);
}

void destroy_clusters(Clusters* clusters){
	//Free all cluster in an array of cluster
	DWORD i;
	for(i=0;i<clusters->size;i++){
		destroy_cluster(clusters->clusters[i]);
	}
	Free(clusters->clusters);
	Free(clusters);
}

void print_cluster(Cluster* cluster,DWORD i){
	printf("    ----Printing cluster %lld----\n",i);
	Object** data=get_objects(cluster);
	DWORD j;
	for(j=0;j<cluster->size;j++){
		printf("%lf,%lf,%lf,%lf\n",data[j]->spatial_coordinates[0],data[j]->spatial_coordinates[1],data[j]->time,data[j]->attribute);
	}
	Free(data);
}

void print_clusters(Clusters* clusters){
	printf("----Printing clusters----\n");
	DWORD i;
	for(i=0;i<clusters->size;i++){
		print_cluster(clusters->clusters[i],i);
	}
}

void print_objects(Object** data,DWORD size){
	DWORD i;
	printf("----Printing objects----\n");
	for(i=0;i<size;i++){
		printf("    ----Printing object %lld----\n",i);
		printf("%lf,%lf,%lf,%lf\n",data[i]->spatial_coordinates[0],data[i]->spatial_coordinates[1],data[i]->time,data[i]->attribute);
	}
}

void remove_from_cluster(Cluster* cluster,Node* previous,Node* current){
	cluster->size-=1;
	//Detach a node current from a link list knowing the previous node.
	if(current==cluster->head){
		cluster->head=current->next;
	}else{
		if(current==cluster->tail){
			previous->next=NULL;
			cluster->tail=previous;
		}else{
			previous->next=current->next;
		}
	}
}

void insert_to_cluster(Cluster* cluster,Node* previous,Node* current){
	cluster->size+=1;
	//Attach the node current after previous node. Maintaining the structure of rest of link list.
	if(previous==NULL){
		current->next=cluster->head;
		cluster->head=current;
	}else{
		if(previous==cluster->tail){
			previous->next=current;
			cluster->tail=current;
		}else{
			current->next=previous->next;
			previous->next=current;
		}
	}
}

void add_to_cluster_front(Cluster* cluster,Object* object){
	if(cluster->size==0){
		//If empty cluster, insert the object as the first object.
		cluster->tail=Calloc(1,Node);
		cluster->tail->object=object;
		cluster->size=1;
		cluster->head=cluster->tail;
		cluster->tail->next=NULL;
		cluster->sum=object->attribute;
	}else{
		//If non-empty cluster, insert the object to the front of cluster (link list)
		Node* head=Calloc(1,Node);
		head->next=cluster->head;
		cluster->head=head;
		cluster->head->object=object;
		cluster->size+=1;
		cluster->sum+=object->attribute;
	}
}

void remove_cluster_front(Cluster* cluster){
	//Remove the front element of the cluster.
	if(cluster->size>0){
		Node* temp=cluster->head;
		cluster->head=temp->next;
		Free(temp);
		cluster->size-=1;
	}
}
