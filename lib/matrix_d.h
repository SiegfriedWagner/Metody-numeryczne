//
// Created by Mateusz Chojnowski on 11.05.2021.
//

#include <stdlib.h>
#include <stdio.h>
#ifndef NUMERYCZNE_MATRIX_D_H
#define NUMERYCZNE_MATRIX_D_H
typedef struct matrix {
    const int rows;
    const int cols;
    double* const * const data;
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
void copy_row(matrix source, matrix target, size_t sourceRow, size_t targetRow);
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
void solve_backward_offset(matrix left, matrix right, matrix output, size_t offset);
matrix solve_equation(matrix A, matrix free_terms);
void invert_matrix_inplace(matrix mat);
void mat_mul_scalar_inplace(matrix mat, double scalar);
double mat_sum(matrix mat);
void mat_neg_inplace(matrix mat);
void mat_sub_mat_h(matrix left, matrix right, matrix output);
double mat_norm(matrix mat);
double vec_norm(matrix vec);
void diagonal(matrix mat, double diagonal_value);
parsing_result from_file(FILE* file);
parsing_code readInt(FILE *file, int *out, char *buffer, unsigned int bufferSize);
parsing_code maybeReadInt(FILE *file, int *out, char *buffer, unsigned int bufferSize);
parsing_code readDouble(FILE *file, double *out, char *buffer, const unsigned int bufferSize);
parsing_code equationFromFile(FILE* file, matrix* mat, matrix* vec);
void print_matrix(matrix mat);
#endif //NUMERYCZNE_MATRIX_D_H
