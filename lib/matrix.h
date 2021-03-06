//
// Created by Mateusz Chojnowski on 30.03.2021.
//

#ifndef MATRIX_H_
#define MATRIX_H_
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

typedef struct matrix {
    const int rows;
    const int cols;
    float* const * const data;
} matrix;

typedef enum parsing_code {
    CORRECT = 0,
    BUFFER_OVERFLOW = 1,
    INVALID_MATRIX_SIZE = 2,
    MATRIX_VALUE_PARSING_ERROR = 3,
    EXCEED_MATRIX_SIZE = 4,
    INVALID_FILE = 5,
    NOT_FOUND = 6
} parsing_code;

typedef struct parsing_result {
    const struct matrix output_matrix;
    const parsing_code output_code;
} parsing_result;

matrix mat_create(int rows, int columns);
matrix mat_create_zero(int rows, int columns);
matrix mat_create_identity(size_t size);
matrix mat_create_copy(matrix mat);
void mat_move(const matrix *from, const matrix *to);
void copy_values(matrix source, matrix target);
void copy_column(matrix source, matrix target, int sourceColumn, int targetColumn);
void destroy_matrix(matrix mat);
matrix mat_mul_mat(matrix a, matrix b);
void mat_mul_mat_h(matrix a, matrix b, matrix output);
matrix transpose(matrix mat);
void transpose_h(matrix transposed, matrix output);
void transpose_inplace(matrix mat);
void swap_cols(matrix mat, int a, int b);
void swap_rows(matrix mat, int a, int b);
void zero_matrix(matrix mat);
// makes LU decomposition - L matrix and U matrix should be created with appropriate size before calling function
void doolitleLU(matrix source, matrix L_out, matrix U_out);
void doolitleLUP(matrix source, matrix L_out, matrix U_out, matrix P_out);
void solve_forward(matrix left, matrix right, matrix output);
void solve_backward(matrix left, matrix right, matrix output);
matrix solve_equation(matrix A, matrix free_terms);
void invert_matrix_inplace(matrix mat);
void mat_mul_scalar_inplace(matrix mat, float scalar);
void mat_neg_inplace(matrix mat);
float mat_norm(matrix mat);
float vec_norm(matrix vec);
void diagonal(matrix mat, float diagonal_value);
parsing_result from_file(FILE* file);
parsing_code readInt(FILE *file, int *out, char *buffer, unsigned int bufferSize);
parsing_code maybeReadInt(FILE *file, int *out, char *buffer, unsigned int bufferSize);
parsing_code readFloat(FILE *file, float *out, char *buffer, unsigned int bufferSize);
parsing_code equationFromFile(FILE* file, matrix* mat, matrix* vec);
void print_matrix(matrix mat);
#endif // MATRIX_H_