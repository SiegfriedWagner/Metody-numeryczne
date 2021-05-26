//
// Created by Mateusz Chojnowski on 9.05.2021.
//
#include <math.h>
#include <assert.h>
#include <float.h>
#include "../lib/matrix_d.h"
#include "../lib/array_d.h"
#include "../lib/interop_d.h"
#include "../lib/svd.h"

int main(int argc, char *argv[]) {
    if (argc != 2)
        printf("Required parameters: matrix file path");
    matrix A = {-1, -1, NULL};
    {
        FILE *file = fopen(argv[1], "r");
        parsing_result result = from_file(file);
        if (result.output_code != CORRECT) {
            printf("Unable to read matrix from file, parsing code: %d", result.output_code);
            return 1;
        }
        mat_move(&result.output_matrix, &A);
        fclose(file);
    }
    print_matrix(A);
    if (A.rows < A.cols) {
        printf("Invalid matrix, rows number should be greater or equal than cols number");
        return 2;
    }
    matrix U = mat_create(A.rows, A.rows);
    matrix Sigma = mat_create(A.rows, A.cols);
    matrix V = mat_create(A.cols, A.cols);
    matrix A_orig = mat_create_copy(A);
    SVD(A, U, V, Sigma);
    matrix Sigma_psudo_inv = mat_create_zero(Sigma.cols, Sigma.rows);
    assert(Sigma.rows >= Sigma.cols);
    for (int i = 0; i < Sigma.cols; ++i) {
        if (fabs(Sigma.data[i][i]) != 0.0)
            Sigma_psudo_inv.data[i][i] = 1.0 / Sigma.data[i][i];
    }
    printf("\nSigma inv\n");
    print_matrix(Sigma_psudo_inv);
    {
        matrix temp = mat_mul_mat(Sigma, Sigma_psudo_inv);
        printf("\nSigma * Sigma_pseudo_inv\n");
        print_matrix(temp);
        destroy_matrix(temp);
    }
    matrix UT = transpose(U);
    matrix temp = mat_mul_mat(V, Sigma_psudo_inv);
    matrix A_puseudo_inv = mat_mul_mat(temp, UT);
    printf("\nA_pseudo_inv\n");
    print_matrix(A_puseudo_inv);
    matrix unit = mat_mul_mat(A_puseudo_inv, A_orig);
    printf("\nA_pseudo_inv * A\n");
    print_matrix(unit);
    destroy_matrix(A_puseudo_inv);
    destroy_matrix(UT);
    destroy_matrix(Sigma_psudo_inv);
    destroy_matrix(U);
    destroy_matrix(Sigma);
    destroy_matrix(V);
    destroy_matrix(A);
    return 0;
}