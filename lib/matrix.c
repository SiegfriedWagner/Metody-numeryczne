//
// Created by Mateusz Chojnowski on 30.03.2021.
//

#include "matrix.h"
const unsigned short bufferSize = 40;

matrix create_matrix(int rows, int columns) {
    // creates matrix with allocated memory but uninitialized
    float **data = malloc(sizeof(float*) * rows);

    data = malloc(sizeof(float*) * rows);
    for (int i = 0; i < rows; ++i) {
        data[i] = malloc(sizeof(float ) * columns);
    }
    matrix mat = { rows, columns, data };
    return mat;
}

matrix create_zero_matrix(int rows, int columns) {
    // creates matrix with allocated memory but uninitialized
    float **data = calloc(rows, sizeof(float*));
    for (int i = 0; i < rows; ++i) {
        data[i] = calloc(columns, sizeof(float));
    }
    matrix mat = {rows, columns, data};
    return mat;
}

matrix copy_matrix(matrix mat) {
    matrix returned = create_matrix(mat.rows, mat.cols);
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            returned.data[row][col] = mat.data[row][col];
        }
    }
    return returned;
}

void copy_values(matrix source, matrix target) {
    assert(source.rows == target.cols);
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

void destroy_matrix(matrix mat) {
    if (mat.data != NULL) {
        for (int i = 0; i < mat.rows; ++i) {
            if (mat.data[i] != NULL) {
                free(mat.data[i]);
            }
        }
        free(mat.data);
    }
    *((int*) &mat.rows) = -1;
    *((int*) &mat.cols) = -1;
    *((float***) &mat.data) = NULL;
}

matrix matmul(matrix a, matrix b) {
    assert(a.cols == b.rows);
    matrix mat = create_matrix(a.rows, b.cols);
    matmul_h(a, b, mat);
    return mat;
}

void matmul_h(matrix a, matrix b, matrix output) {
    assert(a.cols == b.rows && output.rows == a.rows && output.cols == b.cols);
    for (int r = 0; r < output.rows; ++r) {
        for (int c = 0; c < output.cols; ++c) {
            output.data[r][c] = 0.0f;
            for (int i = 0; i < a.cols; ++i) {
                output.data[r][c] += a.data[r][i] * b.data[i][c];
            }
        }
    }
}

matrix transpose(matrix mat) {
    matrix returned = create_matrix(mat.cols, mat.rows);
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            returned.data[col][row] = mat.data[row][col];
        }
    }
    return returned;
}

void transpose_inplace(matrix mat) {
    float temp;
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
    float tmp;
    for (int i = 0; i < mat.rows; ++i) {
        tmp = mat.data[i][a];
        mat.data[i][a] = mat.data[i][b];
        mat.data[i][b] = tmp;
    }
}

void swap_rows(matrix mat, int a, int b) {
    if (a == b)
        return;
    float tmp;
    for (int i = 0; i < mat.cols; ++i) {
        tmp = mat.data[a][i];
        mat.data[a][i] = mat.data[b][i];
        mat.data[b][i] = tmp;
    }
}

void zero_matrix(matrix mat) {
    for (int i = 0; i < mat.rows; ++i)
        for (int j = 0; j < mat.cols; ++j)
            mat.data[i][j] = 0.0f;
}

// makes LU decomposition - L matrix and U matrix should be created with appropriate size before calling function
void doolitleLU(matrix source, matrix L_out, matrix U_out) {
    const int n = source.rows;
    float sum;
    assert(n == L_out.rows && n == U_out.rows);
    assert(n == L_out.cols && n == U_out.cols);
    assert(source.cols == n);
    // fill L matrix with 1 on diagonal
    for (int i = 0; i < source.rows; ++i)
        L_out.data[i][i] = 1.0f;
    for (int k = 0; k < n; ++k) {
        // fill row
        for(int j = k; j < n; ++j) {
            sum = 0.0f;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[k][p] * U_out.data[p][j];
            U_out.data[k][j] = source.data[k][j] - sum;
        }
        // fill column
        for (int i = k + 1; i < n; ++i) {
            sum = 0.0f;
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
    float sum;
    assert(n == L_out.rows && n == U_out.rows);
    assert(n == L_out.cols && n == U_out.cols);
    assert(source.cols == n);
    // fill L matrix with 1 on diagonal
    for (int i = 0; i < source.rows; ++i) {
        L_out.data[i][i] = 1.0f;
        P_out.data[i][i] = 1.0f;
    }
    for (int k = 0; k < n; ++k) {
        // fill row
        for(int j = k; j < n; ++j) {
            sum = 0.0f;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[k][p] * U_out.data[p][j];
            U_out.data[k][j] = source.data[k][j] - sum;
        }
        // find max col index in k row
        float max_value = FLT_MIN;
        int max_index = -1;
        for (int i = k; i < n; ++i) {
            float curr_value = fabs(U_out.data[k][i]);
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
            sum = 0.0f;
            for (int p = 0; p < k; ++p)
                sum += L_out.data[i][p] * U_out.data[p][k];
            L_out.data[i][k] = (source.data[i][k] - sum) / U_out.data[k][k];
        }
    }
}

void solve_forward(matrix left, matrix right, matrix output) {
    assert(right.cols == output.cols && right.rows == output.rows);
    for (int row = 0; row < left.rows; ++row) {
        float b_val = right.data[row][0];
        for (int col = 0; col < row; ++col) {
            b_val -= left.data[row][col] * output.data[col][0];
        }
        output.data[row][0] = b_val / left.data[row][row];
    }
}

void solve_backward(matrix left, matrix right, matrix output) {
    assert(right.cols == output.cols && right.rows == output.rows);
    for (int row = left.rows - 1; row >= 0 ; --row) {
        float y_val = right.data[row][0];
        for (int col = left.cols - 1; col > row; --col) {
            y_val -= left.data[row][col] * output.data[col][0];
        }
        output.data[row][0] = y_val / left.data[row][row];
    }
}

matrix solve_equation(matrix A, matrix free_terms) {
    assert(A.cols == free_terms.rows);
    {
        matrix L = create_zero_matrix(A.rows, A.cols), U = create_zero_matrix(A.rows, A.cols), P = create_zero_matrix(
                A.rows, A.cols);
        doolitleLUP(A, L, U, P);
        matrix vec2 = create_zero_matrix(free_terms.rows, 1);
        solve_forward(L, free_terms, vec2);
        solve_backward(U, vec2, free_terms);
        transpose_inplace(P);
        matmul_h(P, free_terms, vec2);
        destroy_matrix(L);
        destroy_matrix(U);
        destroy_matrix(P);
        return vec2;
    }
}

void invert_matrix_inplace(matrix mat) {
    matrix L = create_zero_matrix(mat.rows, mat.cols), U = create_zero_matrix(mat.rows, mat.cols), P = create_zero_matrix(mat.rows, mat.cols);
    doolitleLUP(mat, L, U, P);
    matrix L_inverted = create_zero_matrix(L.rows, L.cols), U_inverted = create_zero_matrix(U.rows, U.cols);
    matrix supp_vec = create_zero_matrix(L.rows, 1), copy_col_vec = create_zero_matrix(L.rows, 1);
    // invert L and U
    for (int col = 0; col < L.cols; ++col) {
        supp_vec.data[col][0] = 1.0f;
        solve_forward(L, supp_vec, copy_col_vec);
        copy_column(copy_col_vec, L_inverted, 0, col);
        solve_backward(U, supp_vec, copy_col_vec);
        copy_column(copy_col_vec, U_inverted, 0, col);
        supp_vec.data[col][0] = 0.0f;
    }
    transpose_inplace(P);
    matmul_h(P, U_inverted, U);
    matmul_h(U, L_inverted, mat);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(L_inverted);
    destroy_matrix(U_inverted);
    destroy_matrix(supp_vec);
    destroy_matrix(copy_col_vec);
}

float norm(matrix mat) {
    // inf
    float _norm = 0.0f;
    for (int row = 0; row < mat.rows; ++row) {
        float colsum = 0.0f;
        for (int col = 0; col < mat.cols; ++col) {
            colsum += fabs(mat.data[row][col]);
        }
        if (colsum > _norm)
            _norm = colsum;
    }
    return _norm;
}

float vec_norm(matrix vec) {
    assert(vec.cols == 1 || vec.rows == 1);
    float sum = 0.0f;
    for (int row = 0; row < vec.rows; ++row) {
        sum += vec.data[row][0];
    }
    return sum;
}

parsing_result from_file(FILE* file) {
    char buffer[bufferSize];
    unsigned short current_buffer_position = 0;
    int character = fgetc(file);
    while((isalnum(character) || character == '.') && current_buffer_position < bufferSize) {
        buffer[current_buffer_position++] = (char) character;
        character = fgetc(file);
    }
    if(current_buffer_position == bufferSize - 1) {
        parsing_result result = {
                .output_matrix = {-1, -1, NULL},
                .output_code = BUFFER_OVERFLOW
        };
        return result;
    }
    buffer[current_buffer_position++] = '\0';
    int parsed = strtol(buffer, NULL, 10);
    if (parsed <= 0L)
    {
        parsing_result result = {
                .output_matrix = {-1, -1, NULL},
                .output_code = INVALID_MATRIX_SIZE
        };
        return result;
    }
    struct matrix result_matrix = create_matrix(parsed, parsed);
    int row = 0;
    int col = 0;
    current_buffer_position = 0;
    character = fgetc(file);
    while (!feof(file)) {
        while ((isalnum(character) || character == '.') && current_buffer_position < bufferSize) {
            buffer[current_buffer_position++] = (char) character;
            character = fgetc(file);
        }
        if(current_buffer_position == bufferSize - 1) {
            parsing_result result = {
                    .output_matrix = {-1, -1, NULL},
                    .output_code = BUFFER_OVERFLOW
            };
            return result;
        }
        buffer[current_buffer_position] = '\0';
        float parsedf = strtof(buffer, NULL);
        if(parsedf == 0.0f && buffer[0] != '0') { // probably invalid value parsed
            parsing_result result = {
                    .output_matrix = {-1, -1, NULL},
                    .output_code = MATRIX_VALUE_PARSING_ERROR
            };
            return result;
        }
        current_buffer_position = 0;
        if (row >= result_matrix.rows || col >=result_matrix.cols) {
            parsing_result result = {
                    .output_matrix = {-1, -1, NULL},
                    .output_code = EXCEED_MATRIX_SIZE
            };
            return result;
        }
        result_matrix.data[row][col++] = parsedf;
        while(isspace(character)) {
            if(character == '\n') {
                row++;
                col = 0;
            }
            character = getc(file);
        }
    }
    parsing_result result = { .output_matrix = result_matrix, .output_code = CORRECT };
    return result;
}

parsing_code equationFromFile(FILE *file, matrix *mat, matrix *vec)
{
    char buffer[bufferSize];
    int eq_num = 0;
    if (readInt(file, &eq_num, buffer, bufferSize) != CORRECT){
        printf("Error while reading equations num");
        return INVALID_MATRIX_SIZE;
    }
    matrix mat_temp = create_matrix(eq_num, eq_num);
    *((int*) &mat->rows) = mat_temp.rows;
    *((int*) &mat->cols) = mat_temp.cols;
    *((float***) &mat->data) = mat_temp.data;
    matrix vec_temp = create_matrix(eq_num, 1);
    *((int*) &vec->rows) = vec_temp.rows;
    *((int*) &vec->cols) = vec_temp.cols;
    *((float***) &vec->data) = vec_temp.data;
    // fill matrix
    for (int row = 0; row < mat->rows; ++row) {
        for (int col = 0; col < mat->cols; ++col) {
            float val = 0.0f;
            if (readFloat(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                return MATRIX_VALUE_PARSING_ERROR;
            }
            mat->data[row][col] = val;
        }
    }
    // fill b vector
    for (int row = 0; row < vec->rows; ++row) {
        for (int col = 0; col < vec->cols; ++col) {
            float val = 0.0f;
            if (readFloat(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                return MATRIX_VALUE_PARSING_ERROR;
            }
            vec->data[row][col] = val;
        }
    }
    return CORRECT;
}

void printMatrix(matrix mat) {
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            printf("%f ", mat.data[i][j]);
        }
        printf("\n");
    }
}

parsing_code readFloat(FILE* file, float *out, char *buffer, const unsigned  int bufferSize) {
    int character = fgetc(file);
    while(!(isalnum(character) || character == '.' || character == '-' || character == '+')) {
        character = fgetc(file);
    }
    unsigned short current_buffer_position = 0;
    while ((isalnum(character) || character == '.' || character == '-' || character == '+') && current_buffer_position < bufferSize) {
        buffer[current_buffer_position++] = (char) character;
        character = fgetc(file);
    }
    if(current_buffer_position == bufferSize - 1)
        return BUFFER_OVERFLOW;
    buffer[current_buffer_position++] = '\0';
    *out = strtof(buffer, NULL);
    if(*out == 0.0f && buffer[0] != '0') { // probably invalid value parsed
        return MATRIX_VALUE_PARSING_ERROR;
    }
    return CORRECT;
}

parsing_code readInt(FILE *file, int *out, char *buffer, const unsigned int bufferSize) {
    int character = fgetc(file);
    while(!(isalnum(character) || character == '.' || character == '-' || character == '+'))
        character = fgetc(file);
    unsigned short current_buffer_position = 0;
    while ((isalnum(character) || character == '.' || character == '-' || character == '+') && current_buffer_position < bufferSize) {
        buffer[current_buffer_position++] = (char) character;
        character = fgetc(file);
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
