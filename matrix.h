/*
 * Copyright (C) 2016, Northwestern University.
 * This file defines a matrix data structure.
 * Matrix operations are implemented.
 * See also matrix_functions.c for implementations.
*/

#ifndef KRIG_MATRIX_H
#define KRIG_MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "clustertype.h"

/*
 * Definition for matrix in 2 dimensional array pointer form.
*/
typedef struct{
	DWORD n_row;
	DWORD n_column;
	DTYPE **matrix;
}Matrix;
/*
 * Initialize a matrix by giving its column and row numbers.
 * n_row: Number of rows.
 * n_column: Number of columns.
 * Return: A matrix structure with memory space allocated.
*/
extern Matrix* create_matrix(DWORD n_row,DWORD n_column);
/*
 * Destroy a matrix. Free all memory space.
 * matrix: The pointer to the matrix you want to destroy.
*/
extern void destroy_matrix(Matrix* matrix);
/*
 * Add two matrix together.
 * m1: The first matrix.
 * m2: The second matrix.
 * Return: A new matrix that is the sum of the two matrices.
*/
extern Matrix* matrix_addition(Matrix* m1,Matrix* m2);
/*
 * Multiply two matrix together. Automatically determine if naive or strassen algorithm is used based on matrix size.
 * m1: The first matrix.
 * m2: The second matrix.
 * Return: A new matrix that is the product of the two matrices.
*/
extern Matrix* matrix_multiplication(Matrix* m1,Matrix* m2);
/*
 * Invert an upper triangular matrix. (Do not check input)
 * m: The matrix to be inverted.
 * Return: A new matrix that is the inversion of m.
*/
extern Matrix* upper_triangle_inversion(Matrix* m);
/*
 * Cholesky decomposition for a positive definite matrix and output its inversion.
 * m: A positive definite matrix to be inverted.
 * Return: A new matrix that is the inversion of m.
*/
extern Matrix* cholesky_decomposition(Matrix* m);
/*
 * Transpose a matrix.
 * m: The matrix to be transposed.
 * Return: A new matrix that is the transpose of m.
*/
extern Matrix* matrix_transpose(Matrix* m);
/*
 * Find lower upper permutation of a matrix m.
 * m: The matrix to be decomposed.
 * Return: 3 matrices. The first one is L, the second one is U, and the third one is P. We have Pm=LU.
*/
extern Matrix** lower_upper_permutation(Matrix* m);
/*
 * Solve a linear system of equation using lower upper permutation.
 * Gamma: The coefficients of linear system on LHS.
 * gamma: The result of linear system on RHS.
 * Return: The value of variables at LHS. In the end Gamma*variables=gamma.
*/
extern Matrix* solve_linear_system(Matrix* Gamma,Matrix* gamma);
/*
 * Multiply two matrix together using naive matrix multiplication algorithm.
 * m1: The first matrix.
 * m2: The second matrix.
 * Return: A new matrix that is the product of the two matrices.
*/
extern Matrix* naive_multiplication(Matrix* m1,Matrix* m2);
/*
 * Multiply two matrix together using strassen matrix multiplication algorithm.
 * A: The first matrix.
 * B: The second matrix.
 * Return: A new matrix that is the product of the two matrices.
*/
extern Matrix* strassen_multiplication(Matrix* A,Matrix* B);
extern void resize_matrix(Matrix* m,DWORD n_row,DWORD n_column);
extern void print_matrix(Matrix* m);
extern Matrix* vector_self_multiplication(Matrix* m);
extern void add_to_matrix(Matrix* m1,Matrix* m2);

#endif
