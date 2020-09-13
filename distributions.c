/*
 * Copyright (C) 2016, Northwestern University.
 * This file contains functions for generating random data order or random data from statistical distrubtion.
 * See also random.h
*/


#include "cluster.h"
#include "random.h"

DWORD poisson(DTYPE lambda){
	DWORD k=0;
	DTYPE L=-lambda;
	DTYPE p=0;
	do{
		k++;
		//Logrithm is used to avoid underflow.
		p+=log(genrand_real2());
	} while(p>L);
	return --k;
}

unsigned int* random_ints(int n,int size){
	unsigned int* indices=Calloc(n,unsigned int);
	unsigned int i,temp;
	//Create an array of size n for positive integer sequence 1,...,n
	for(i=0;i<(unsigned int)n;i++){
		indices[i]=i;
	}
	//Shuffle the array
	for(i=(unsigned int)n-1;i>0;i--){
		unsigned long r=genrand_int32()%i;
		temp=indices[i];
		indices[i]=indices[r];
		indices[r]=temp;
	}
	//Keep the first size number of integers.
	indices=(unsigned int*)realloc(indices,size*sizeof(unsigned int));
	return indices;
}

void shuffle_objects(Object** data,DWORD size){
	Object* temp;
	DWORD i;
	for(i=size-1;i>0;i--){
		DWORD r=genrand_int32()%(i+1);
		temp=data[i];
		data[i]=data[r];
		data[r]=temp;
	}
}
