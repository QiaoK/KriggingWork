#include "matrix.h"

/*
 *===========================Description===========================
 * Initialize a matrix given row and column sizes.
 *===========================Parameters============================
 * n_row: Number of rows.
 * n_column: Number of columns
 *===========================Return================================
 * The pointer to the matrix created.
*/
Matrix* create_matrix(DWORD n_row,DWORD n_column){
	Matrix* result=Calloc(1,Matrix);
	result->n_row=n_row;
	result->n_column=n_column;
	DWORD i=0;
	result->matrix=Calloc(n_row,DTYPE*);
	for(i=0;i<n_row;i++){
		result->matrix[i]=Calloc(n_column,DTYPE);
	}
	return result;
}
/*
 * Truncate matrix into chuncks.
*/
void resize_matrix(Matrix* m,DWORD n_row,DWORD n_column){
	DWORD i;
	for(i=0;i<n_row;i++){
		m->matrix[i]=(DTYPE*)realloc(m->matrix[i],sizeof(DTYPE)*n_column);
	}
	m->matrix=(DTYPE**)realloc(m->matrix,sizeof(DTYPE*)*n_row);
}

/*
 * Free all pointers for a matrix.
*/

void destroy_matrix(Matrix* matrix){
	if(matrix==NULL){
		return;
	}
	DWORD i=0;
	for(i=0;i<matrix->n_row;i++){
		Free(matrix->matrix[i]);
	}
	Free(matrix);
}
/*
 * Print the matrix into command line.
*/
void print_matrix(Matrix* m){
	if(m==NULL){
		printf("Empty matrix\n");
		return;
	}
	DWORD i,j;
	for(i=0;i<m->n_row;i++){
		for(j=0;j<m->n_column;j++){
			printf("%lf ",m->matrix[i][j]);
		}
		printf("\n");
	}
}
/*
 * Add two matrices together, return a new matrix.
*/
Matrix* matrix_addition(Matrix* m1,Matrix* m2){
	if(m1->n_column!=m2->n_column||m1->n_row!=m2->n_row){
		return NULL;
	}
	Matrix* result=create_matrix(m1->n_row,m1->n_column);
	DWORD i,j;
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m1->n_column;j++){
			result->matrix[i][j]=m1->matrix[i][j]+m2->matrix[i][j];
		}
	}
	return result;
}

/*
 * Add m2 to matrix m1.
*/

void add_to_matrix(Matrix* m1,Matrix* m2){
	if(m1->n_column!=m2->n_column||m1->n_row!=m2->n_row){
		return;
	}
	DWORD i,j;
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m1->n_column;j++){
			m1->matrix[i][j]+=m2->matrix[i][j];
		}
	}
}

Matrix* matrix_subtraction(Matrix* m1,Matrix* m2){
	if(m1->n_column!=m2->n_column||m1->n_row!=m2->n_row){
		return NULL;
	}
	Matrix* result=create_matrix(m1->n_row,m1->n_column);
	DWORD i,j;
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m1->n_column;j++){
			result->matrix[i][j]=m1->matrix[i][j]-m2->matrix[i][j];
		}
	}
	return result;
}

Matrix* read_upper_triangle_solution(Matrix* Gamma,Matrix* gamma){
	DWORD i,j;
	Matrix* result=create_matrix(gamma->n_row,1);
	for(i=gamma->n_row-1;i>=0;i--){
		result->matrix[i][0]=gamma->matrix[i][0];
		for(j=gamma->n_row-1;j>i;j--){
			result->matrix[i][0]-=result->matrix[j][0]*Gamma->matrix[i][j];
		}
		if(Gamma->matrix[i][i]!=0){
			result->matrix[i][0]/=Gamma->matrix[i][i];
		}else{
			destroy_matrix(result);
			return NULL;
		}
	}
	return result;
}

Matrix* read_lower_triangle_solution(Matrix* Gamma,Matrix* gamma){
	DWORD i,j;
	Matrix* result=create_matrix(gamma->n_row,1);
	for(i=0;i<gamma->n_row;i++){
		result->matrix[i][0]=gamma->matrix[i][0];
		for(j=0;j<i;j++){
			result->matrix[i][0]-=result->matrix[j][0]*Gamma->matrix[i][j];
		}
		if(Gamma->matrix[i][i]!=0){
			result->matrix[i][0]/=Gamma->matrix[i][i];
		}else{
			destroy_matrix(result);
			return NULL;
		}
	}
	return result;
}
/*
Matrix* solve_linear_system(Matrix* Gamma,Matrix* gamma){
	Matrix *m1,*m2,*m3;
	Matrix** LUP=lower_upper_permutation(Gamma);
	Matrix* left=matrix_multiplication(LUP[2],gamma);
	destroy_matrix(LUP[2]);
	m1=matrix_transpose(LUP[0]);
	destroy_matrix(LUP[0]);
	m3=upper_triangle_inversion(m1);
	destroy_matrix(m1);
	m2=matrix_transpose(m3);
	destroy_matrix(m3);
	m1=matrix_multiplication(m2,left);
	destroy_matrix(m2);
	destroy_matrix(left);
	m3=upper_triangle_inversion(LUP[1]);
	destroy_matrix(LUP[1]);
	m2=matrix_multiplication(m3,m1);
	destroy_matrix(m1);
	destroy_matrix(m3);
	Free(LUP);
	return m2;
}
*/

Matrix* solve_linear_system(Matrix* Gamma,Matrix* gamma){
	Matrix *m1,*m2;
	Matrix** LUP=lower_upper_permutation(Gamma);
	Matrix* left=Calloc(1,Matrix);
	left->n_row=gamma->n_row;
	left->n_column=1;
	left->matrix=Calloc(left->n_row,DTYPE*);
	DWORD i;
	for(i=0;i<gamma->n_row;i++){
		left->matrix[i]=gamma->matrix[(DWORD)LUP[3]->matrix[i][0]];
	}
	destroy_matrix(LUP[2]);
	destroy_matrix(LUP[3]);
	m1=read_lower_triangle_solution(LUP[0],left);
	destroy_matrix(LUP[0]);
	Free(left->matrix);
	Free(left);
	if(m1!=NULL){
		m2=read_upper_triangle_solution(LUP[1],m1);
	}else{
		m2=NULL;
	}
	destroy_matrix(LUP[1]);
	destroy_matrix(m1);
	Free(LUP);
	return m2;
}

Matrix* expand_to_power(Matrix* m,DWORD n){
	DTYPE exponent=log(n)/log(2);
	if(floorf(exponent)!=exponent){
		n=(DWORD)pow(2,floorf(exponent)+1);
	}
	//printf("n=%ld\n",n);
	Matrix* result=create_matrix(n,n);
	DWORD i,j;
	for(i=0;i<m->n_row;i++){
		memcpy(result->matrix[i],m->matrix[i],m->n_column*sizeof(DTYPE));
		for(j=m->n_column;j<n;j++){
			result->matrix[i][j]=0;
		}
	}
	for(i=m->n_row;i<n;i++){
		for(j=0;j<n;j++){
			result->matrix[i][j]=0;
		}
	}
	return result;
}

Matrix* resize(Matrix* m,DWORD n_row,DWORD n_column){
	Matrix* result=create_matrix(n_row,n_column);
	DWORD i;
	for(i=0;i<result->n_row;i++){
		memcpy(result->matrix[i],m->matrix[i],n_column*sizeof(DTYPE));
	}
	return result;
}

Matrix* naive_multiplication(Matrix* m1,Matrix* m2){
	if(m1->n_column!=m2->n_row){
		return NULL;
	}
	Matrix* result=create_matrix(m1->n_row,m2->n_column);
	DWORD i,j,k;
	for(i=0;i<m1->n_row;i++){
		for(j=0;j<m2->n_column;j++){
			result->matrix[i][j]=0;
			for(k=0;k<m1->n_column;k++){
				result->matrix[i][j]+=m1->matrix[i][k]*m2->matrix[k][j];
			}
		}
	}
	return result;
}

Matrix* sub_matrix(Matrix* m,DWORD i,DWORD j){
	Matrix* result=Calloc(1,Matrix);
	DWORD n=m->n_row/2;
	result->n_row=n;
	result->n_column=n;
	result->matrix=Calloc(n,DTYPE*);
	DWORD row_jump=(i-1)*n;
	DWORD column_jump=(j-1)*n;
	DWORD k;
	for(k=0;k<n;k++){
		result->matrix[k]=m->matrix[k+row_jump]+column_jump;
	}
	return result;
}

Matrix* strassen_multiplication(Matrix* A,Matrix* B){
	if(A->n_row<=64){
		return naive_multiplication(A,B);
	}
	Matrix* A11=sub_matrix(A,1,1);
	Matrix* A12=sub_matrix(A,1,2);
	Matrix* A21=sub_matrix(A,2,1);
	Matrix* A22=sub_matrix(A,2,2);
	Matrix* B11=sub_matrix(B,1,1);
	Matrix* B12=sub_matrix(B,1,2);
	Matrix* B21=sub_matrix(B,2,1);
	Matrix* B22=sub_matrix(B,2,2);
	//Compute M1
	Matrix* temp1=matrix_addition(A11,A22);
	Matrix* temp2=matrix_addition(B11,B22);
	Matrix* M1=matrix_multiplication(temp1,temp2);
	destroy_matrix(temp1);
	destroy_matrix(temp2);
	//Compute M2
	temp1=matrix_addition(A21,A22);
	Matrix* M2=matrix_multiplication(temp1,B11);
	destroy_matrix(temp1);
	//Compute M3
	temp1=matrix_subtraction(B12,B22);
	Matrix* M3=matrix_multiplication(A11,temp1);
	destroy_matrix(temp1);
	//Compute M4
	temp1=matrix_subtraction(B21,B11);
	Matrix* M4=matrix_multiplication(A22,temp1);
	destroy_matrix(temp1);
	//Compute M5
	temp1=matrix_addition(A11,A12);
	Matrix* M5=matrix_multiplication(temp1,B22);
	//Compute M6
	temp1=matrix_subtraction(A21,A11);
	temp2=matrix_addition(B11,B12);
	Matrix* M6=matrix_multiplication(temp1,temp2);
	destroy_matrix(temp1);
	destroy_matrix(temp2);
	//Compute M7
	temp1=matrix_subtraction(A12,A22);
	temp2=matrix_addition(B21,B22);
	Matrix* M7=matrix_multiplication(temp1,temp2);
	destroy_matrix(temp1);
	destroy_matrix(temp2);
	//Compute C
	Matrix* C=create_matrix(A->n_row,A->n_column);
	DWORD i,j;
	DWORD n=A->n_row/2;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			C->matrix[i][j]=M1->matrix[i][j]+M4->matrix[i][j]-M5->matrix[i][j]+M7->matrix[i][j];
			C->matrix[i][j+n]=M3->matrix[i][j]+M5->matrix[i][j];
			C->matrix[i+n][j]=M2->matrix[i][j]+M4->matrix[i][j];
			C->matrix[i+n][j+n]=M1->matrix[i][j]-M2->matrix[i][j]+M3->matrix[i][j]+M6->matrix[i][j];
		}
	}
	destroy_matrix(M1);
	destroy_matrix(M2);
	destroy_matrix(M3);
	destroy_matrix(M4);
	destroy_matrix(M5);
	destroy_matrix(M6);
	destroy_matrix(M7);
	Free(A11);
	Free(A12);
	Free(A21);
	Free(A22);
	Free(B11);
	Free(B12);
	Free(B21);
	Free(B22);
	return C;
}

Matrix* vector_self_multiplication(Matrix* m){
	if(m->n_column!=1){
		return NULL;
	}
	Matrix* result=create_matrix(m->n_row,m->n_row);
	DWORD i,j;
	for(i=0;i<m->n_row;i++){
		for(j=0;j<m->n_row;j++){
			result->matrix[i][j]=m->matrix[i][0]*m->matrix[j][0];
		}
	}
	return result;
}

Matrix* matrix_multiplication(Matrix* m1,Matrix* m2){
	if(m1->n_column!=m2->n_row){
		return NULL;
	}
	DWORD n;
	if(m1->n_row>=m1->n_column){
		if(m1->n_row>=m2->n_column){
			n=m1->n_row;
		}else{
			n=m2->n_column;
		}
	}else{
		if(m1->n_column>=m2->n_column){
			n=m1->n_column;
		}else{
			n=m2->n_column;
		}
	}
	Matrix* A=expand_to_power(m1,n);
	Matrix* B=expand_to_power(m2,n);
	Matrix* C;
	if(A->n_row>64){
		C=strassen_multiplication(A,B);
	}else{
		C=naive_multiplication(A,B);
	}
	destroy_matrix(A);
	destroy_matrix(B);
	Matrix* result=resize(C,m1->n_row,m2->n_column);
	destroy_matrix(C);
	return result;
}

Matrix* upper_triangle_inversion(Matrix* m){
	if(m->n_column!=m->n_row){
		return NULL;
	}
	Matrix* result=create_matrix(m->n_row,m->n_column);
	DWORD i,j,k;
	DTYPE multiplier;
	for(i=0;i<m->n_row;i++){
		if(m->matrix[i][i]==0){
			return NULL;
		}
		for(j=0;j<m->n_column;j++){
			result->matrix[i][j]=0;
		}
	}
	for(i=0;i<m->n_row&&m->matrix[i][i]!=0;i++){
		result->matrix[i][i]=1;
		for(k=0;k<=i;k++){
			result->matrix[k][i]/=m->matrix[i][i];
		}
		for(j=i+1;j<m->n_column;j++){
			multiplier=m->matrix[i][j];
			for(k=0;k<=i;k++){
				result->matrix[k][j]-=multiplier*result->matrix[k][i];
			}
		}
	}
	return result;
}

Matrix* cholesky_decomposition(Matrix* m){
	if(m->n_column!=m->n_row){
		return NULL;
	}
	Matrix* result=create_matrix(m->n_row,m->n_column);
	DWORD i,j,k;
	for(j=m->n_column-1;j>=0;j--){
		for(i=m->n_row-1;i>j;i--){
			result->matrix[i][j]=0;
		}
		result->matrix[j][j]=m->matrix[j][j];
		for(k=m->n_column-1;k>j;k--){
		//printf("%lf\n",result->matrix[j][k]);
			result->matrix[j][j]-=result->matrix[j][k]*result->matrix[j][k];
		}
		//printf("%lf\n",m->matrix[j][j]);
		result->matrix[j][j]=sqrt(result->matrix[j][j]);
		for(i=j-1;i>=0;i--){
			result->matrix[i][j]=m->matrix[i][j];
			for(k=m->n_column-1;k>j;k--){
				result->matrix[i][j]-=result->matrix[i][k]*result->matrix[j][k];
			}
			if(result->matrix[j][j]==0){
				result->matrix[i][j]=0;
			}else{
				result->matrix[i][j]/=result->matrix[j][j];
			}
		}
	}
	return result;
}

void swap_row(Matrix* m,DWORD x,DWORD y){
	DTYPE* temp;
	temp=m->matrix[x];
	m->matrix[x]=m->matrix[y];
	m->matrix[y]=temp;
}

Matrix** lower_upper_permutation(Matrix* m){
	if(m->n_column!=m->n_row){
		return NULL;
	}
	Matrix** result=Calloc(4,Matrix*);
	Matrix* P=create_matrix(m->n_row,m->n_column);
	Matrix* U=create_matrix(m->n_row,m->n_column);
	Matrix* L=create_matrix(m->n_row,m->n_column);
	Matrix* CP=create_matrix(m->n_row,1);
	result[0]=L;
	result[1]=U;
	result[2]=P;
	result[3]=CP;
	DWORD i,j,k;
	DTYPE multiplier,temp;
	for(i=0;i<m->n_row;i++){
		for(j=0;j<m->n_column;j++){
			U->matrix[i][j]=m->matrix[i][j];
			P->matrix[i][j]=0;
			L->matrix[i][j]=0;
		}
		L->matrix[i][i]=1;
		P->matrix[i][i]=1;
		CP->matrix[i][0]=i;
	}
	for(i=0;i<m->n_column;i++){
		j=i;
		while(j<m->n_row&&U->matrix[j][i]==0){
			j++;
		}
		if(j==m->n_row){
			break;
		}
		if(j!=i){
			swap_row(P,i,j);
			swap_row(U,i,j);
			swap_row(CP,i,j);
			for(k=0;k<i;k++){
				temp=L->matrix[i][k];
				L->matrix[i][k]=L->matrix[j][k];
				L->matrix[j][k]=temp;
			}
		}
		for(j=i+1;j<m->n_row;j++){
			multiplier=U->matrix[j][i]/U->matrix[i][i];
			L->matrix[j][i]=multiplier;
			for(k=0;k<m->n_column;k++){
				U->matrix[j][k]-=multiplier*U->matrix[i][k];
			}
		}
	}
	return result;
}

Matrix* matrix_transpose(Matrix* m){
	Matrix* result=create_matrix(m->n_column,m->n_row);
	DWORD i,j;
	for(i=0;i<result->n_row;i++){
		for(j=0;j<result->n_column;j++){
			result->matrix[i][j]=m->matrix[j][i];
		}
	}
	return result;
}
