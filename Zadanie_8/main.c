//
// Created by Mateusz Chojnowski on 20.04.2021.
//
#include <assert.h>
#include <math.h>
#include "../lib/matrix.h"
#include "../lib/array.h"


void QRdecomposition(matrix A, matrix Q, matrix R) {
    assert(A.rows == Q.rows && A.rows == R.rows);
    assert(A.cols == Q.cols && A.cols == R.cols);
    assert(A.rows == A.cols);
    matrix A_k = mat_create_copy(A);
    matrix calc_helping_matrix = mat_create(A.rows, A.cols);
    diagonal(Q, 1.0f);
    array x_k = create_zero_array(A.rows);
    matrix H = mat_create_identity(A.rows);
    for (int k = 0; k < A.rows - 1; ++k) {
        for (int i = 0; i < A.rows; ++i) {
            x_k.data[i] = 0.0f;
        }
        diagonal(H, 1.0f);
        float alpha = 0.0f;
        for (int i = 0; i < A.rows - k; ++i) {
            x_k.data[i] = A_k.data[i + k][k];
            alpha += x_k.data[i] * x_k.data[i];
        }
        alpha = sqrtf(alpha);
        if (A.data[k][k] > 0)
            alpha = -alpha;
        // Based on wikipedia
        // x_k == v
        x_k.data[0] -= alpha; // u = x - alpha * e_1
        // v = u / ||u||
        float dotprod2 = 2 / arrdot(x_k, x_k);
        for (int i = k; i < A.rows; ++i) {
            for (int j = k; j < A.rows; ++j) {
                // H = 1 - 2 * v * v
                H.data[i][j] = H.data[i][j] - (x_k.data[i - k] * dotprod2 * x_k.data[j - k]);
            }
        }
        copy_values(Q, calc_helping_matrix);
        mat_mul_mat_h(calc_helping_matrix, H, Q);
        copy_values(A_k, calc_helping_matrix);
        mat_mul_mat_h(H, calc_helping_matrix, A_k);
    }
    transpose_h(Q, calc_helping_matrix); // Q^T
    mat_mul_mat_h(calc_helping_matrix, A, R);
    destroy_matrix(H);
    destroy_array(x_k);
    destroy_matrix(calc_helping_matrix);
    destroy_matrix(A_k);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Missing parameter: matrix file path\n");
        printf("Optional parameter: maximum difference between QR and RQ (per matrix value)\n");
        return 1;
    }
    float maximumError = 0.0f;
    if (argc == 3)
        maximumError = strtof(argv[2], NULL);
    FILE *file;
    file = fopen(argv[1], "r");
    parsing_result parsingResult = from_file(file);
    fclose(file);
    if (parsingResult.output_code != CORRECT) {
        printf("Error during parsing, error code: %d", parsingResult.output_code);
        return 2;
    }
    matrix A = parsingResult.output_matrix;
    matrix Q = mat_create(A.rows, A.cols);
    matrix R = mat_create(A.rows, A.cols);
    matrix QR = mat_create(A.rows, A.cols);
    int step = 0;
    int iterate = 1;
    while(iterate) {
        printf("\nStep: %d\n", step);
        zero_matrix(Q);
        zero_matrix(R);
        printf("\nA\n");
        print_matrix(A);
        QRdecomposition(A, Q, R);
        printf("\nQ\n");
        print_matrix(Q);
        printf("\nR\n");
        print_matrix(R);
        mat_mul_mat_h(Q, R, QR);
        printf("\nQR\n");
        print_matrix(QR);
        mat_mul_mat_h(R, Q, A); // RQ
        step++;
        iterate = 1;
        // check diagonal equality end condition
        for (int i = 0; i < A.rows; ++i) {
            if (A.data[i][i] == QR.data[i][i]){
                iterate = 0;
            }
            else {
                iterate = 1;
                break;
            }
        }
        // check maximum error end condition
        if(maximumError > 0 && iterate == 1) {
            for (int i = 0; i < A.rows; ++i) {
                for (int j = 0; j < A.cols; ++j) {
                    if(fabs(A.data[i][j] - QR.data[i][j]) < maximumError) {
                        iterate = 0;
                    }
                    else
                    {
                        iterate = 1;
                        break;
                    }
                }
                if(iterate == 1) break;
            }
        }
    }
    printf("\nEigenvalues\n");
    for (int i = 0; i < A.rows; ++i) {
        printf("%f ", A.data[i][i]);
    }

    destroy_matrix(QR);
    destroy_matrix(R);
    destroy_matrix(Q);
    destroy_matrix(A);
    return 0;
}