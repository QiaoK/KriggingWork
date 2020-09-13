/*
 * This file contains basic functions for random processes.
 * See distrubtion.c
 * Use implmentation of MT19937.
*/

#ifndef KRIG_RANDOM_H
#define KRIG_RANDOM_H

#include "cluster.h"
/*
 * A C-program for MT19937, with initialization improved 2002/1/26.
 *  Coded by Takuji Nishimura and Makoto Matsumoto.
 * Check liscence in random.c
*/
extern void init_genrand(unsigned long s);
extern void init_by_array(unsigned long init_key[], int key_length);
extern unsigned long genrand_int32(void);
extern double genrand_real2(void);

/*
 * Generator for poisson Distribution.
 * Use MT19937.
 * lambda: Poisson mean.
 * RETURN: Random value in discrete Poisson distribution.
*/
extern DWORD poisson(DTYPE lambda);

/*
 * Fisher-Yates shuffle method for objects
 * Use MT19937.
 * data: Array of objects to be shuffled.
 * size: Size of object array.
*/

extern void shuffle_objects(Object** data,DWORD size);
/*
 * Fisher-Yates shuffle method to random select without replacement.
 * Use MT19937.
 * n: The upper bound value.
 * size: The size of array for result.
 * Return: Size number of integers in range 0 to n-1 without repetition. 
*/
extern unsigned int* random_ints(int n,int size);

#endif
