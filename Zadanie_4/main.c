//
// Created by Mateusz Chojnowski on 17.03.2021.
//

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

const unsigned short bufferSize = 40;

typedef struct matrix {
    const int rows;
    const int cols;
    float** data;
} matrix;


matrix create_matrix(int rows, int columns) {
    // creates matrix with allocated memory but uninitialized
    matrix mat = {rows, columns, NULL};
    mat.data = malloc(sizeof(float*) * rows);
    for (int i = 0; i < rows; ++i) {
        mat.data[i] = malloc(sizeof(float ) * columns);
    }
    return mat;
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

int main(int argc, char **argv) {
    FILE* file;
    if(argc != 2) {
        printf("Missing argument: input file ");
        return 1;
    }

    file = fopen(argv[1], "r");
    parsing_result result = fromFile(file);
    if(result.output_code == CORRECT) {
        printf("Matrix:\n");
        printMatrix(result.output_matrix);
        matrix L = create_matrix(result.output_matrix.rows, result.output_matrix.cols);
        matrix U = create_matrix(result.output_matrix.rows, result.output_matrix.cols);
        matrix P = create_matrix(result.output_matrix.rows, result.output_matrix.cols);
        zero_matrix(L);
        zero_matrix(U);
        zero_matrix(P);
        doolitleLUP(result.output_matrix, L, U, P);
        //doolitleLU(result.output_matrix, L, U);
        printf("L:\n");
        printMatrix(L);
        printf("U':\n");
        printMatrix(U);
        printf("P:\n");
        printMatrix(P);
        matrix recreated = matmul(L, U);
        printf("L*U'");
        printMatrix(recreated);
        matrix LUP = matmul(recreated, P);
        printf("L*U'*P\n");
        printMatrix(LUP);
        destroy_matrix(LUP);
        destroy_matrix(L);
        destroy_matrix(U);
    }
    fclose(file);
    destroy_matrix(result.output_matrix);
}