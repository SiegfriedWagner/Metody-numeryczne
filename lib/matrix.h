//
// Created by mateu on 30.03.2021.
//

#ifndef MATRIX_H_
#define MATRIX_H_
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>

typedef struct matrix {
    const int rows;
    const int cols;
    float** data;
} matrix;

typedef enum parsing_code {
    CORRECT = 0,
    BUFFER_OVERFLOW = 1,
    INVALID_MATRIX_SIZE = 2,
    MATRIX_VALUE_PARSING_ERROR = 3,
    EXCEED_MATRIX_SIZE = 4
} parsing_code;

typedef struct parsing_result {
    const struct matrix output_matrix;
    const parsing_code output_code;
} parsing_result;

matrix create_matrix(int rows, int columns);
matrix create_zero_matrix(int rows, int columns);
matrix copy_matrix(matrix mat);
//void scale_equation(matrix )
void destroy_matrix(matrix mat);
matrix matmul(matrix a, matrix b);
void transposeInplace(matrix mat);
void swap_cols(matrix mat, int a, int b);
void swap_rows(matrix mat, int a, int b);
void zero_matrix(matrix mat);
// makes LU decomposition - L matrix and U matrix should be created with appropriate size before calling function
void doolitleLU(matrix source, matrix L_out, matrix U_out);
void doolitleLUP(matrix source, matrix L_out, matrix U_out, matrix P_out);
void solve_forward(matrix left, matrix right, matrix output);
void solve_backward(matrix left, matrix right, matrix output);
parsing_result fromFile(FILE* file);
parsing_code readInt(FILE *file, int *out, char *buffer, unsigned int bufferSize);
parsing_code readFloat(FILE *file, float *out, char *buffer, unsigned int bufferSize);
void printMatrix(matrix mat);
#endif // MATRIX_H_