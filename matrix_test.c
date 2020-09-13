/*
*  Copyright (C) 2016, Northwestern University.
*/
#include <stdio.h>
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#define EPS 0.000001

/*
 * Unit tests for matrix implementation.
 * All test cases are examples for how to use the matrix implementation.
 * Matrix operation is the low level part Kriging interpolation.
*/

BOOLEAN test_identity(Matrix* m){
	DWORD i,j;
	for(i=0;i<m->n_row;i++){
		for(j=0;j<m->n_column;j++){
			if(i==j){
				if(fabs(m->matrix[i][j]-1)>EPS){
					return FALSE;
				}
			}else{
				if(fabs(m->matrix[i][j])>EPS){
					return FALSE;
				}
			}
		}
	}
	return TRUE;
}

BOOLEAN test_equality(Matrix* m1,Matrix* m2){
	if(m1->n_row!=m2->n_row||m1->n_column!=m2->n_column){
		return FALSE;
	}
	DWORD i,j;
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m1->n_column;j++){
			if(fabs(m1->matrix[i][j]-m2->matrix[i][j])>EPS){
				return FALSE;
			}
		}
	}
	return TRUE;
}


int main(void){
	Matrix *m1,*m2,*result,*m3,*gamma;
	Matrix** LUP;
	DWORD i,j;

	// Test for matrix multiplication.
	m1=create_matrix(2,2);
	m1->matrix[0][0]=1;
	m1->matrix[1][0]=3;
	m1->matrix[0][1]=2;
	m1->matrix[1][1]=4;
	result=matrix_multiplication(m1,m1);
	if(result->matrix[0][0]!=7.){
		printf("Test 1 failed: matrix[0][0]=%lf!=7\n",result->matrix[0][0]);
		return -1;
	}
	if(result->matrix[0][1]!=10.){
		printf("Test 1 failed: matrix[0][1]=%lf!=10\n",result->matrix[0][1]);
		return -1;
	}
	if(result->matrix[1][0]!=15.){
		printf("Test 1 failed: matrix[1][0]=%lf!=15\n",result->matrix[1][0]);
		return -1;
	}
	if(result->matrix[1][1]!=22.){
		printf("Test 1 failed: matrix[1][1]=%lf!=22\n",result->matrix[1][1]);
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(result);
	
	// Test for upper triangle matrix inversion.
	m1=create_matrix(3,3);
	m1->matrix[0][0]=1;
	m1->matrix[0][1]=1;
	m1->matrix[0][2]=1;
	m1->matrix[1][0]=0;
	m1->matrix[1][1]=1;
	m1->matrix[1][2]=1;
	m1->matrix[2][0]=0;
	m1->matrix[2][1]=0;
	m1->matrix[2][2]=1;
	result=upper_triangle_inversion(m1);
	m2=matrix_multiplication(m1,result);
	//print_matrix(result);
	if(!test_identity(m2)){
		printf("Test 2 failed: Inversion not correct.");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(result);

	// Test for upper triangle matrix inversion.
	m1=create_matrix(3,3);
	m1->matrix[0][0]=2;
	m1->matrix[0][1]=9;
	m1->matrix[0][2]=11;
	m1->matrix[1][0]=0;
	m1->matrix[1][1]=4;
	m1->matrix[1][2]=12;
	m1->matrix[2][0]=0;
	m1->matrix[2][1]=0;
	m1->matrix[2][2]=15;
	//print_matrix(m1);
	m2=upper_triangle_inversion(m1);
	//print_matrix(m2);
	result=matrix_multiplication(m1,m2);
	//print_matrix(result);
	if(!test_identity(result)){
		printf("Test 3 failed: Inversion not correct.");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(result);

	// Test for symmetric matrix inversion with cholesky decomposition.
	m1=create_matrix(3,3);
	m1->matrix[0][0]=4;
	m1->matrix[0][1]=12;
	m1->matrix[0][2]=-16;
	m1->matrix[1][0]=12;
	m1->matrix[1][1]=37;
	m1->matrix[1][2]=-43;
	m1->matrix[2][0]=-16;
	m1->matrix[2][1]=-43;
	m1->matrix[2][2]=98;
	//print_matrix(m1);
	result=cholesky_decomposition(m1);
	//print_matrix(result);
	m2=upper_triangle_inversion(result);
	//print_matrix(matrix_multiplication(result,m2));
	m3=matrix_transpose(m2);
	destroy_matrix(result);
	result=matrix_multiplication(m3,m2);
	//print_matrix(result);
	destroy_matrix(m2);
	m2=matrix_multiplication(m1,result);
	//print_matrix(m2);
	if(!test_identity(m2)){
		printf("Test 4 failed: Inversion not correct.");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(m3);
	destroy_matrix(result);

	//Test for LUR decomposition without row swapping.
	m1=create_matrix(3,3);
	m1->matrix[0][0]=0;
	m1->matrix[0][1]=-7;
	m1->matrix[0][2]=0;
	m1->matrix[1][0]=-3;
	m1->matrix[1][1]=0;
	m1->matrix[1][2]=6;
	m1->matrix[2][0]=5;
	m1->matrix[2][1]=-1;
	m1->matrix[2][2]=0;
	LUP=lower_upper_permutation(m1);
	m2=matrix_multiplication(LUP[0],LUP[1]);
	m3=matrix_multiplication(LUP[2],m1);
	if(!test_equality(m2,m3)){
		printf("Test 5 failed: LUR decomposition does not add up to PA=LU");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(m3);
	destroy_matrix(LUP[0]);
	destroy_matrix(LUP[1]);
	destroy_matrix(LUP[2]);
	destroy_matrix(LUP[3]);
	Free(LUP);

	//Test for LUR decomposition with row swapping.
	m1=create_matrix(3,3);
	m1->matrix[0][0]=10;
	m1->matrix[0][1]=-7;
	m1->matrix[0][2]=0;
	m1->matrix[1][0]=10;
	m1->matrix[1][1]=-7;
	m1->matrix[1][2]=2;
	m1->matrix[2][0]=20;
	m1->matrix[2][1]=-7;
	m1->matrix[2][2]=5;
	LUP=lower_upper_permutation(m1);
	m2=matrix_multiplication(LUP[0],LUP[1]);
	m3=matrix_multiplication(LUP[2],m1);
	if(!test_equality(m2,m3)){
		printf("Test 6 failed: LUR decomposition does not add up to PA=LU");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(m3);
	destroy_matrix(LUP[0]);
	destroy_matrix(LUP[1]);
	destroy_matrix(LUP[2]);
	destroy_matrix(LUP[3]);
	Free(LUP);

	//Test for LUR decomposition without row swapping, a 4 by 4 matrix with.
	m1=create_matrix(4,4);
	m1->matrix[0][0]=10;
	m1->matrix[0][1]=-7;
	m1->matrix[0][2]=2;
	m1->matrix[0][3]=3;
	m1->matrix[1][0]=10;
	m1->matrix[1][1]=-7;
	m1->matrix[1][2]=1;
	m1->matrix[1][3]=4;
	m1->matrix[2][0]=12;
	m1->matrix[2][1]=-7;
	m1->matrix[2][2]=5;
	m1->matrix[2][3]=5;
	m1->matrix[3][0]=10;
	m1->matrix[3][1]=-7;
	m1->matrix[3][2]=2;
	m1->matrix[3][3]=12;
	LUP=lower_upper_permutation(m1);
	m2=matrix_multiplication(LUP[0],LUP[1]);
	m3=matrix_multiplication(LUP[2],m1);
	if(!test_equality(m2,m3)){
		printf("Test 7 failed: LUR decomposition does not add up to PA=LU");
		return -1;
	}
	//print_matrix(LUP[2]);
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(m3);
	destroy_matrix(LUP[0]);
	destroy_matrix(LUP[1]);
	destroy_matrix(LUP[2]);
	destroy_matrix(LUP[3]);
	Free(LUP);

	//Test for large matrix multiplication
	m1=create_matrix(129,330);
	m2=create_matrix(330,175);
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m1->n_column;j++){
			m1->matrix[i][j]=(i+j)*(i-j);
		}
	}
	for(i=0;i<m2->n_row;i++){
		for(j=0;j<m2->n_column;j++){
			m2->matrix[i][j]=(i+j+1)*(i-j+1);
		}
	}
	result=naive_multiplication(m1,m2);
	m3=matrix_multiplication(m1,m2);
	if(!test_equality(result,m3)){
		printf("Test 8 failed: Multiplication results are not equal.");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(m3);
	destroy_matrix(result);

	//Test for very large matrix multiplication
	m1=create_matrix(523,1156);
	m2=create_matrix(1156,367);
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m1->n_column;j++){
			m1->matrix[i][j]=(i+j)*(i-j);
		}
	}
	for(i=0;i<m2->n_row;i++){
		for(j=0;j<m2->n_column;j++){
			m2->matrix[i][j]=(i+j+1)*(i-j+1);
		}
	}
	result=naive_multiplication(m1,m2);
	m3=matrix_multiplication(m1,m2);
	if(!test_equality(result,m3)){
		printf("Test 9 failed: Multiplication results are not equal.");
		return -1;
	}
	destroy_matrix(m1);
	destroy_matrix(m2);
	destroy_matrix(m3);
	destroy_matrix(result);

	//Test for solving linear system, not automated.
	m1=create_matrix(3,3);
	m1->matrix[0][0]=0;
	m1->matrix[0][1]=0.4751;
	m1->matrix[0][2]=1;
	m1->matrix[1][0]=0.4751;
	m1->matrix[1][1]=0;
	m1->matrix[1][2]=1;
	m1->matrix[2][0]=1;
	m1->matrix[2][1]=1;
	m1->matrix[2][2]=0;
	gamma=create_matrix(3,1);
	gamma->matrix[0][0]=0.4629;
	gamma->matrix[1][0]=0.2256;
	gamma->matrix[2][0]=1;
	m2=solve_linear_system(m1,gamma);
	//printf("check.\n");
	//print_matrix(m2);
	destroy_matrix(m1);
	destroy_matrix(m2);
	//print_matrix(LUP[2]);
	//Test for singular matrix

	m1=create_matrix(3,3);
	m1->matrix[0][0]=0;
	m1->matrix[0][1]=2;
	m1->matrix[0][2]=3;
	m1->matrix[1][0]=1;
	m1->matrix[1][1]=0;
	m1->matrix[1][2]=3;
	m1->matrix[2][0]=1;
	m1->matrix[2][1]=2;
	m1->matrix[2][2]=0;

	gamma=create_matrix(3,1);
	gamma->matrix[0][0]=10;
	gamma->matrix[1][0]=3;
	gamma->matrix[2][0]=5;

	LUP=lower_upper_permutation(m1);

	/*print_matrix(LUP[0]);
	print_matrix(LUP[1]);
	print_matrix(LUP[2]);
	print_matrix(LUP[3]);*/
	m2=matrix_multiplication(LUP[0],LUP[1]);
	m3=matrix_multiplication(LUP[2],m1);
	if(!test_equality(m2,m3)){
		printf("Test 10 failed: LUR decomposition does not add up to PA=LU");
		return -1;
	}
	destroy_matrix(m2);
	//print_matrix(m1);
	m2=solve_linear_system(m1,gamma);
	
	destroy_matrix(m1);
	destroy_matrix(gamma);

	//print_matrix(m2);
	destroy_matrix(m2);
	destroy_matrix(LUP[0]);
	destroy_matrix(LUP[1]);
	destroy_matrix(LUP[2]);
	destroy_matrix(LUP[3]);
	Free(LUP);

	printf("Test finished.\n");
	return 0;
}
