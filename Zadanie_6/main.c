//
// Created by Mateusz Chojnowski on 31.03.2021.
//

#include "../lib/matrix.h"
#include <math.h>

void scaleMatrixWithResultsInplace(matrix mat, matrix result) {
    // scales matrix by finding largest element in mat row and dividing every element in row by max and dividing
    // result row by max
    assert(mat.rows == mat.cols);
    assert(mat.cols == result.rows);
    assert(result.cols == 1);
    for (int row = 0; row < mat.rows; ++row) {
        float max = 0.0f;
        for (int col = 0; col < mat.cols; ++col) {
            max = fmaxf(max, mat.data[row][col]);
        }
        for (int col = 0; col < mat.cols; ++col) {
            mat.data[row][col] /= max;
        }
        result.data[row][0] /= max;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Missing argument: input file ");
        return 1;
    }
    FILE *file = fopen(argv[1], "r");
    matrix mat = { -1, -1, NULL};
    matrix b = { -1, -1, NULL };
    if (equationFromFile(file, &mat, &b) != CORRECT) {
        printf("Error during file parsing\n");
        return 2;
    }
    printf("A\n");
    print_matrix(mat);
    matrix A = mat_create_copy(mat);
    printf("\nB\n");
    print_matrix(b);
    // calcualtions
    matrix L = mat_create_zero(mat.rows, mat.cols), U = mat_create_zero(mat.rows, mat.cols), P = mat_create_zero(
            mat.rows, mat.cols);
    doolitleLUP(mat, L, U, P);
    matrix L_inverted = mat_create_zero(L.rows, L.cols), U_inverted = mat_create_zero(U.rows, U.cols);
    matrix supp_vec1 = mat_create_zero(L.rows, 1), supp_vec2 = mat_create_zero(L.rows, 1);
    printf("\nL\n");
    print_matrix(L);
    printf("\nU\n");
    print_matrix(U);
    printf("\nP'\n");
    print_matrix(P);
    transpose_inplace(P);
    // invert L and U
    for (int col = 0; col < L.cols; ++col) {
        supp_vec1.data[col][0] = 1.0f;
        solve_forward(L, supp_vec1, supp_vec2);
        copy_column(supp_vec2, L_inverted, 0, col);
        solve_backward(U, supp_vec1, supp_vec2);
        copy_column(supp_vec2, U_inverted, 0, col);
        supp_vec1.data[col][0] = 0.0f;
    }
    printf("\nU^(-1)\n");
    print_matrix(U_inverted);
    printf("\nL^(-1)\n");
    print_matrix(L_inverted);
    transpose_inplace(P);
    // solve equation
    solve_forward(L, b, supp_vec1);
    solve_backward(U, supp_vec1, supp_vec2); // supp_vec2 = x_prim
    mat_mul_mat_h(P, supp_vec2, supp_vec1); // supp_vec1 = x
    // calculate A^(-1)
    mat_mul_mat_h(P, U_inverted, U);
    mat_mul_mat_h(U, L_inverted, mat); // mat = A^(-1)
    // print result
    printf("\nA^(-1)\n");
    print_matrix(mat);
    matrix test = mat_mul_mat(A, mat);
    printf("\nA * A^(-1)\n");
    print_matrix(test);
    printf("\ncond(A) = %f\n", mat_norm(A) * mat_norm(mat));
    printf("\ncond_estimated(A) = %f\n", mat_norm(A) * vec_norm(supp_vec1) / vec_norm(b));
    mat_mul_mat_h(A, supp_vec1, supp_vec2); // A*x = b
    printf("\nA*x\n");
    print_matrix(supp_vec2);
    // scaled
    copy_values(A, mat);
    scaleMatrixWithResultsInplace(mat, b);
    matrix A_scaled = mat_create_copy(mat);
    printf("\nA (scaled)\n");
    print_matrix(mat);
    printf("\nB (scaled)\n");
    print_matrix(b);
    // calcualtions
    zero_matrix(L);
    zero_matrix(U);
    zero_matrix(P);
    doolitleLUP(mat, L, U, P);
    zero_matrix(L_inverted);
    zero_matrix(U_inverted);
    zero_matrix(supp_vec1);
    zero_matrix(supp_vec2);
    // invert L and U
    for (int col = 0; col < L.cols; ++col) {
        supp_vec1.data[col][0] = 1.0f;
        solve_forward(L, supp_vec1, supp_vec2);
        copy_column(supp_vec2, L_inverted, 0, col);
        solve_backward(U, supp_vec1, supp_vec2);
        copy_column(supp_vec2, U_inverted, 0, col);
        supp_vec1.data[col][0] = 0.0f;
    }
    transpose_inplace(P);
    // solve equation
    solve_forward(L, b, supp_vec1);
    solve_backward(U, supp_vec1, supp_vec2); // supp_vec2 = x_prim
    mat_mul_mat_h(P, supp_vec2, supp_vec1); // supp_vec1 = x
    // calculate A^(-1)
    mat_mul_mat_h(P, U_inverted, U);
    mat_mul_mat_h(U, L_inverted, mat); // mat = A^(-1)
    printf("\nA^(-1)\n");
    print_matrix(mat);
    mat_mul_mat_h(A_scaled, mat, test);
    printf("\nA * A^(-1)\n");
    print_matrix(test);
    printf("\ncond(A) = %f\n", mat_norm(A_scaled) * mat_norm(mat));
    printf("\ncond_estimated(A) = %f\n", mat_norm(A_scaled) * vec_norm(supp_vec1) / vec_norm(b));
    mat_mul_mat_h(A_scaled, supp_vec1, supp_vec2); // A*x = b
    printf("\nA*x\n");
    print_matrix(supp_vec2);

    destroy_matrix(test);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(L_inverted);
    destroy_matrix(U_inverted);
    destroy_matrix(supp_vec1);
    destroy_matrix(supp_vec2);
}