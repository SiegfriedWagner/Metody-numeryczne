//
// Created by Mateusz Chojnowski on 30.03.2021.
//

#include "matrix.h"

const unsigned short bufferSize = 40;


matrix create_matrix(int rows, int columns) {
    // creates matrix with allocated memory but uninitialized
    matrix mat = {rows, columns, NULL};
    mat.data = malloc(sizeof(float*) * rows);
    for (int i = 0; i < rows; ++i) {
        mat.data[i] = malloc(sizeof(float ) * columns);
    }
    return mat;
}

matrix create_zero_matrix(int rows, int columns) {
    // creates matrix with allocated memory but uninitialized
    matrix mat = {rows, columns, NULL};
    mat.data = calloc(rows, sizeof(float*));
    for (int i = 0; i < rows; ++i) {
        mat.data[i] = calloc(columns, sizeof(float));
    }
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

void destroy_matrix(matrix mat) {
    if (mat.data != NULL) {
        for (int i = 0; i < mat.rows; ++i) {
            if (mat.data[i] != NULL) {
                free(mat.data[i]);
                mat.data[i] = NULL;
            }
        }
        free(mat.data);
    }
    mat.data = NULL;
}

matrix matmul(matrix a, matrix b) {
    assert(a.cols == b.rows);
    matrix mat = create_matrix(a.rows, b.cols);
    for (int r = 0; r < mat.rows; ++r) {
        for (int c = 0; c < mat.cols; ++c) {
            mat.data[r][c] = 0.0f;
            for (int i = 0; i < a.cols; ++i) {
                mat.data[r][c] += a.data[r][i] * b.data[i][c];
            }
        }
    }
    return mat;
}

void transposeInplace(matrix mat) {
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

parsing_result fromFile(FILE* file) {
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