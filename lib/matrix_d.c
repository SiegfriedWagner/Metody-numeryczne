//
// Created by Mateusz Chojnowski on 11.05.2021.
//
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "matrix_d.h"

const unsigned short bufferSize = 40;
static int last_read_character = '\0';

matrix mat_create(int rows, int columns) {
    assert(rows > 0);
    assert(columns > 0);
    // creates matrix with allocated memory but uninitialized
    double **data = malloc(sizeof(double*) * rows);
    for (int i = 0; i < rows; ++i) {
        data[i] = malloc(sizeof(double) * columns);
    }
    matrix mat = { rows, columns, data };
    return mat;
}

matrix mat_create_zero(int rows, int columns) {
    assert(rows > 0);
    assert(columns > 0);
    // creates matrix with allocated memory but uninitialized
    double **data = calloc(rows, sizeof(double*));
    for (int i = 0; i < rows; ++i) {
        data[i] = calloc(columns, sizeof(double));
    }
    matrix mat = {rows, columns, data};
    return mat;
}

matrix mat_create_identity(size_t size) {
    matrix returned = mat_create_zero(size, size);
    for (int i = 0; i < size; ++i) {
        returned.data[i][i] = 1.0;
    }
    return returned;
}

matrix mat_create_copy(matrix mat) {
    matrix returned = mat_create(mat.rows, mat.cols);
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            returned.data[row][col] = mat.data[row][col];
        }
    }
    return returned;
}

void mat_move(const matrix *from, const matrix *to) {
    assert(from->cols > 0 && from->rows > 0);
    destroy_matrix(*to);
    *(int*) &to->rows = (int) from->rows;
    *(int*) &to->cols = (int) from->cols;
    *(double***) &(to->data) = (double**) from->data;
    *(double***) &(from->data) = NULL;
}

void copy_values(matrix source, matrix target) {
    assert(source.rows == target.rows);
    assert(source.cols == target.cols);
    for (int i = 0; i < source.rows; ++i) {
        for (int j = 0; j < source.cols; ++j) {
            target.data[i][j] = source.data[i][j];
        }
    }
}

void copy_column(matrix source, matrix target, int sourceColumn, int targetColumn) {
    assert(sourceColumn >= 0 && targetColumn >= 0);
    assert(sourceColumn < source.cols && targetColumn < target.cols);
    assert(source.rows == target.rows);
    for (int row = 0; row < source.rows; ++row) {
        target.data[row][targetColumn] = source.data[row][sourceColumn];
    }
}

void copy_row(matrix source, matrix target, size_t sourceRow, size_t targetRow) {
    assert(sourceRow >= 0 && targetRow >= 0);
    assert(sourceRow < source.rows && targetRow < target.rows);
    assert(source.cols == target.cols);
    for (int col = 0; col < source.rows; ++col)
        target.data[targetRow][col] = source.data[sourceRow][col];
}

void destroy_matrix(matrix mat) {
    if (mat.data != NULL) {
        for (int i = 0; i < mat.rows; ++i) {
            if (mat.data[i] != NULL) {
                free(mat.data[i]);
            }
        }
        free((double **) mat.data);
    }
    *((int*) &mat.rows) = -1;
    *((int*) &mat.cols) = -1;
    *((double***) &mat.data) = NULL;
}

matrix mat_mul_mat(matrix a, matrix b) {
    assert(a.cols == b.rows);
    matrix mat = mat_create(a.rows, b.cols);
    mat_mul_mat_h(a, b, mat);
    return mat;
}

void mat_mul_mat_h(matrix a, matrix b, matrix output) {
    // assert(a.cols == b.rows && output.rows == a.rows && output.cols == b.cols);
    double sum = 0;
    for (int r = 0; r < output.rows; ++r) {
        for (int c = 0; c < output.cols; ++c) {
            sum = 0.0;
            for (int i = 0; i < a.cols; ++i) {
                double temp = a.data[r][i] * b.data[i][c];
                sum += temp;
            }
            output.data[r][c] = sum;
        }
    }
}

void mat_sub_mat_h(matrix left, matrix right, matrix output) {
    assert(left.cols == right.cols && left.cols == output.cols);
    assert(left.rows == right.rows && left.cols == output.cols);
    for(int i = 0; i < left.rows; ++i)
        for (int j = 0; j < left.cols; ++j)
            output.data[i][j] = left.data[i][j] - right.data[i][j];
}

matrix transpose(matrix mat) {
    matrix returned = mat_create(mat.cols, mat.rows);
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            returned.data[col][row] = mat.data[row][col];
        }
    }
    return returned;
}

void transpose_h(matrix transposed, matrix output) {
    assert(transposed.cols == output.rows);
    assert(transposed.rows == output.cols);
    for (int row = 0; row < transposed.rows; ++row) {
        for (int col = 0; col < transposed.cols; ++col) {
            output.data[col][row] = transposed.data[row][col];
        }
    }
}

void transpose_inplace(matrix mat) {
    double temp;
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < row; ++col) {
            temp = mat.data[row][col];
            mat.data[row][col] = mat.data[col][row];
            mat.data[col][row] = temp;
        }
    }
}

void swap_cols(matrix mat, int a, int b) {
    if (a == b)
        return;
    double tmp;
    for (int i = 0; i < mat.rows; ++i) {
        tmp = mat.data[i][a];
        mat.data[i][a] = mat.data[i][b];
        mat.data[i][b] = tmp;
    }
}

void swap_rows(matrix mat, int a, int b) {
    if (a == b)
        return;
    double tmp;
    for (int i = 0; i < mat.cols; ++i) {
        tmp = mat.data[a][i];
        mat.data[a][i] = mat.data[b][i];
        mat.data[b][i] = tmp;
    }
}

void zero_matrix(matrix mat) {
    for (int i = 0; i < mat.rows; ++i)
        for (int j = 0; j < mat.cols; ++j)
            mat.data[i][j] = 0.0;
}

// makes LU decomposition - L matrix and U matrix should be created with appropriate size before calling function
void doolitleLU(matrix source, matrix L_out, matrix U_out) {
    const int n = source.rows;
    double sum;
    assert(n == L_out.rows && n == U_out.rows);
    assert(n == L_out.cols && n == U_out.cols);
    assert(source.cols == n);
    // fill L matrix with 1 on diagonal
    for (int i = 0; i < source.rows; ++i)
        L_out.data[i][i] = 1.0;
    for (int k = 0; k < n; ++k) {
        // fill row
        for(int j = k; j < n; ++j) {
            sum = 0.0;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[k][p] * U_out.data[p][j];
            U_out.data[k][j] = source.data[k][j] - sum;
        }
        // fill column
        for (int i = k + 1; i < n; ++i) {
            sum = 0.0;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[i][p] * U_out.data[p][k];
            L_out.data[i][k] = (source.data[i][k] - sum) / U_out.data[k][k];
        }
    }
}

void doolitleLUP(matrix source, matrix L_out, matrix U_out, matrix P_out) {
    // makes LUP decomposition - L, P and U matrix should be created with appropriate size and zeroed before calling function
    // WARNING: changes source matrix
    const int n = source.rows;
    double sum;
    assert(n == L_out.rows && n == U_out.rows);
    assert(n == L_out.cols && n == U_out.cols);
    assert(source.cols == n);
    // fill L matrix with 1 on diagonal
    for (int i = 0; i < source.rows; ++i) {
        L_out.data[i][i] = 1.0;
        P_out.data[i][i] = 1.0;
    }
    for (int k = 0; k < n; ++k) {
        // fill row
        for(int j = k; j < n; ++j) {
            sum = 0.0;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[k][p] * U_out.data[p][j];
            U_out.data[k][j] = source.data[k][j] - sum;
        }
        // find max col index in k row
        double max_value = DBL_MIN;
        int max_index = -1;
        for (int i = k; i < n; ++i) {
            double curr_value = fabs(U_out.data[k][i]);
            if (curr_value > max_value) {
                max_value = curr_value;
                max_index = i;
            }
        }
        if(max_index > 0) {
            swap_cols(source, k, max_index);
            swap_cols(U_out, k, max_index);
            swap_rows(P_out, k, max_index);
        }
        // fill column
        for (int i = k + 1; i < n; ++i) {
            sum = 0.0;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[i][p] * U_out.data[p][k];
            L_out.data[i][k] = (source.data[i][k] - sum) / U_out.data[k][k];
        }
    }
}

void solve_forward(matrix left, matrix right, matrix output) {
    assert(right.cols == output.cols && right.rows == output.rows);
    for (int row = 0; row < left.rows; ++row) {
        double b_val = right.data[row][0];
        for (int col = 0; col < row; ++col) {
            b_val -= left.data[row][col] * output.data[col][0];
        }
        output.data[row][0] = b_val / left.data[row][row];
    }
}

void solve_backward(matrix left, matrix right, matrix output) {
    solve_backward_offset(left, right, output, 0);
}

void solve_backward_offset(matrix left, matrix right, matrix output, size_t offset) {
    assert(right.cols == output.cols && right.rows == output.rows);
    for (int row = left.rows - 1 - offset; row >= 0 ; --row) {
        double y_val = right.data[row][0];
        for (int col = left.cols - 1; col > row; --col) {
            y_val -= left.data[row][col] * output.data[col][0];
        }
        output.data[row][0] = y_val / left.data[row][row];
    }
}

matrix solve_equation(matrix A, matrix free_terms) {
    assert(A.cols == free_terms.rows);
    {
        matrix L = mat_create_zero(A.rows, A.cols), U = mat_create_zero(A.rows, A.cols), P = mat_create_zero(
                A.rows, A.cols);
        doolitleLUP(A, L, U, P);
        matrix vec2 = mat_create_zero(free_terms.rows, 1);
        solve_forward(L, free_terms, vec2);
        solve_backward(U, vec2, free_terms);
        transpose_inplace(P);
        mat_mul_mat_h(P, free_terms, vec2);
        destroy_matrix(L);
        destroy_matrix(U);
        destroy_matrix(P);
        return vec2;
    }
}

void invert_matrix_inplace(matrix mat) {
    matrix L = mat_create_zero(mat.rows, mat.cols), U = mat_create_zero(mat.rows, mat.cols), P = mat_create_zero(
            mat.rows, mat.cols);
    doolitleLUP(mat, L, U, P);
    matrix L_inverted = mat_create_zero(L.rows, L.cols), U_inverted = mat_create_zero(U.rows, U.cols);
    matrix supp_vec = mat_create_zero(L.rows, 1), copy_col_vec = mat_create_zero(L.rows, 1);
    // invert L and U
    for (int col = 0; col < L.cols; ++col) {
        supp_vec.data[col][0] = 1.0;
        solve_forward(L, supp_vec, copy_col_vec);
        copy_column(copy_col_vec, L_inverted, 0, col);
        solve_backward(U, supp_vec, copy_col_vec);
        copy_column(copy_col_vec, U_inverted, 0, col);
        supp_vec.data[col][0] = 0.0;
    }
    transpose_inplace(P);
    mat_mul_mat_h(P, U_inverted, U);
    mat_mul_mat_h(U, L_inverted, mat);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(L_inverted);
    destroy_matrix(U_inverted);
    destroy_matrix(supp_vec);
    destroy_matrix(copy_col_vec);
}

void mat_mul_scalar_inplace(matrix mat, double scalar) {
    for (int i = 0; i < mat.rows; ++i)
        for (int j = 0; j < mat.cols; ++j)
            mat.data[i][j] *= scalar;
}

double mat_sum(matrix mat) {
    double sum = 0.0;
    for (int i = 0; i < mat.rows; ++i)
        for (int j = 0; j < mat.cols; ++j)
            sum += fabs(mat.data[i][j]);
    return sum;
}

void mat_neg_inplace(matrix mat) {
    for (int i = 0; i < mat.rows; ++i)
        for (int j = 0; j < mat.cols; ++j)
            mat.data[i][j] = -mat.data[i][j];
}

double mat_norm(matrix mat) {
    // 1
    double _norm = 0.0;
    for (int row = 0; row < mat.rows; ++row) {
        double colsum = 0.0;
        for (int col = 0; col < mat.cols; ++col) {
            colsum += fabs(mat.data[row][col]);
        }
        if (colsum > _norm)
            _norm = colsum;
    }
    return _norm;
}

double vec_norm(matrix vec) {
    assert(vec.cols == 1);
    double sum = 0.0;
    for (int row = 0; row < vec.rows; ++row) {
        sum += fabs(vec.data[row][0]);
    }
    return sum;
}

void diagonal(matrix mat, double diagonal_value) {
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            if (row != col)
                mat.data[row][col] = 0.0;
            else
                mat.data[row][col] = diagonal_value;
        }
    }
}

parsing_result invalid_result(parsing_code code) {
    parsing_result result = {
            .output_matrix = {-1, -1, NULL},
            .output_code = code
    };
    return result;
}

parsing_result from_file(FILE* file) {
    char buffer[bufferSize];
    int parsed_rows = 0;
    int parsed_cols = 0;
    if (file == NULL)
        return invalid_result(INVALID_FILE);
    parsing_code code = readInt(file, &parsed_rows, buffer, bufferSize);
    if(code != CORRECT)
        return invalid_result(code);
    if (parsed_rows <= 0L)
        return invalid_result(INVALID_MATRIX_SIZE);
    parsed_cols = parsed_rows;
    code = maybeReadInt(file, &parsed_cols, buffer, bufferSize);
    if (code != CORRECT && code != NOT_FOUND)
        return invalid_result(code);
    if (parsed_cols <= 0L)
        return invalid_result(INVALID_MATRIX_SIZE);
    matrix result_matrix = mat_create(parsed_rows, parsed_cols);
    int row = 0;
    int col = 0;
    double parsedf = 0.0;
    while (!feof(file)) {
        code = readDouble(file, &parsedf, buffer, bufferSize);
        if (code != CORRECT)
        {
            break;
        }
        if (row >= result_matrix.rows || col >=result_matrix.cols) {
            code = EXCEED_MATRIX_SIZE;
            break;
        }
        result_matrix.data[row][col++] = parsedf;
        while(isspace(last_read_character)) {
            if(last_read_character == '\n') {
                row++;
                col = 0;
            }
            last_read_character = getc(file);
        }
    }
    if (code != CORRECT)
    {
        destroy_matrix(result_matrix);
        parsing_result result = { .output_matrix = {-1, -1, NULL}, .output_code = code };
        return result;
    }
    else
    {
        parsing_result result = { .output_matrix = result_matrix, .output_code = CORRECT };
        return result;
    }
}

parsing_code equationFromFile(FILE *file, matrix *mat, matrix *vec) {
    char buffer[bufferSize];
    int eq_num = 0;
    if (readInt(file, &eq_num, buffer, bufferSize) != CORRECT){
        printf("Error while reading equations num");
        return INVALID_MATRIX_SIZE;
    }
    matrix mat_temp = mat_create(eq_num, eq_num);
    matrix vec_temp = mat_create(eq_num, 1);
    // fill matrix
    for (int row = 0; row < mat_temp.rows; ++row) {
        for (int col = 0; col < mat_temp.cols; ++col) {
            double val = 0.0;
            if (readDouble(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                destroy_matrix(mat_temp);
                destroy_matrix(vec_temp);
                return MATRIX_VALUE_PARSING_ERROR;
            }
            mat_temp.data[row][col] = val;
        }
    }
    // fill b vector
    for (int row = 0; row < vec_temp.rows; ++row) {
        for (int col = 0; col < vec_temp.cols; ++col) {
            double val = 0.0;
            if (readDouble(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                destroy_matrix(mat_temp);
                destroy_matrix(vec_temp);
                return MATRIX_VALUE_PARSING_ERROR;
            }
            vec_temp.data[row][col] = val;
        }
    }
    mat_move(&mat_temp, mat);
    mat_move(&vec_temp, vec);
    return CORRECT;
}

void print_matrix(matrix mat) {
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            printf("%e ", mat.data[i][j]);
        }
        printf("\n");
    }
}

parsing_code readDouble(FILE* file, double *out, char *buffer, const unsigned  int bufferSize) {
    while(!(isalnum(last_read_character) || last_read_character == '.' || last_read_character == '-' || last_read_character == '+')) {
        last_read_character = fgetc(file);
    }
    unsigned short current_buffer_position = 0;
    while ((isalnum(last_read_character) || last_read_character == '.' || last_read_character == '-' || last_read_character == '+') && current_buffer_position < bufferSize) {
        buffer[current_buffer_position++] = (char) last_read_character;
        last_read_character = fgetc(file);
    }
    if(current_buffer_position == bufferSize - 1)
        return BUFFER_OVERFLOW;
    buffer[current_buffer_position++] = '\0';
    *out = strtod(buffer, NULL);
    if(*out == 0.0 && buffer[0] != '0') { // probably invalid value parsed
        return MATRIX_VALUE_PARSING_ERROR;
    }
    return CORRECT;
}

parsing_code readInt(FILE *file, int *out, char *buffer, const unsigned int bufferSize) {
    while(!(isalnum(last_read_character) || last_read_character == '-' || last_read_character == '+'))
        last_read_character = fgetc(file);
    unsigned short current_buffer_position = 0;
    while ((isalnum(last_read_character) || last_read_character == '-' || last_read_character == '+') && current_buffer_position < bufferSize) {
        buffer[current_buffer_position++] = (char) last_read_character;
        last_read_character = fgetc(file);
    }
    if(current_buffer_position == bufferSize - 1)
        return BUFFER_OVERFLOW;
    buffer[current_buffer_position++] = '\0';
    *out = strtol(buffer, NULL, 10);
    if(*out == 0 && buffer[0] != '0') { // probably invalid value parsed
        return MATRIX_VALUE_PARSING_ERROR;
    }
    return CORRECT;
}

parsing_code maybeReadInt(FILE *file, int *out, char *buffer, const unsigned int bufferSize) {
    while(!(isalnum(last_read_character) || last_read_character == '+' || last_read_character == '-' || last_read_character == '\n'))
        last_read_character = fgetc(file);
    if (last_read_character == '\n')
        return NOT_FOUND;
    return readInt(file, out, buffer, bufferSize);
}
