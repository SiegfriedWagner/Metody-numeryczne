//
// Created by mateu on 21.04.2021.
//
#include "../lib/matrix.h"
#include "../lib/array.h"

void matmul_scalar_inplace(matrix mat, float scalar) {
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            mat.data[row][col] *= scalar;
        }
    }
}


int main(int argc, char *argv[]) {
    if(argc < 3) {
        printf("Required parameters: matrix file path and eigenvalues (separated with space)\n");
        return 1;
    }
    FILE *file = fopen(argv[1], "r");
    parsing_result parsingResult = from_file(file);
    fclose(file);
    if(parsingResult.output_code != CORRECT) {
        printf("Error during parsing, error code: %d", parsingResult.output_code);
        return 2;
    }
    matrix A = parsingResult.output_matrix;
    matrix zero_vec = mat_create_zero(A.cols, 1);
    array eigenvalues = create_array(argc - 2);
    for (int i = 2; i < argc; ++i)
        eigenvalues.data[i - 2] = strtof(argv[i], NULL);
    printf("A\n");
    print_matrix(A);
    printf("\neigenvalues\n");
    print_array(eigenvalues);
    // algorytm 2
    matrix L = mat_create(A.rows, A.cols), U = mat_create(A.rows, A.cols), P = mat_create(A.rows, A.cols), A_lambda = mat_create(
            A.rows, A.cols);
    matrix x_prim = mat_create(A.cols, 1);
    for (int i = 0; i < eigenvalues.size; ++i) {
        float eigenvalue = eigenvalues.data[i];
        printf("\nEigenvalue: %f", eigenvalue);
        zero_matrix(L);
        zero_matrix(U);
        zero_matrix(P);
        copy_values(A, A_lambda);
        zero_matrix(zero_vec);
        for (int j = 0; j < A_lambda.rows; ++j) {
            A_lambda.data[j][j] -= eigenvalue; // - epislon);
        }
        doolitleLUP(A_lambda, L, U, P);
        x_prim.data[x_prim.rows - 1][0] = 1.0f;
        for (int row = U.rows - 2; row >= 0 ; --row) {
            float y_val = 0;
            for (int col = U.cols - 1; col > row; --col) {
                y_val -= U.data[row][col] * x_prim.data[col][0];
            }
            x_prim.data[row][0] = y_val / U.data[row][row];
        }
        transpose_inplace(P);
        mat_mul_mat_h(P, x_prim, zero_vec);
        printf("\nEigenvector (x_i)\n");
        print_matrix(zero_vec);
        mat_mul_mat_h(A, zero_vec, x_prim);
        printf("\nA * x_i\n");
        print_matrix(x_prim);
        matmul_scalar_inplace(zero_vec, eigenvalue);
        printf("\neigenvalue * x_i\n");
        print_matrix(zero_vec);
    }
    destroy_matrix(x_prim);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(A_lambda);
    destroy_matrix(zero_vec);
    destroy_array(eigenvalues);
    destroy_matrix(A);
}