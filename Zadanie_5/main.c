//
// Created by Mateusz Chojnowski on 24.03.2021.
//

#include "../lib/matrix.h"

int main(int argc, char *argv[]) {
    if(argc != 2) {
        printf("Missing argument: input file ");
        return 1;
    }
    FILE *file = fopen(argv[1], "r");
    matrix mat = {-1, -1, NULL };
    matrix vec1 = { -1, -1, NULL };
    if (equationFromFile(file, &mat, &vec1) != CORRECT) {
        printf("Error during parsing input file");
        return 2;
    }
    matrix L = create_zero_matrix(mat.rows, mat.cols), U = create_zero_matrix(mat.rows, mat.cols), P = create_zero_matrix(mat.rows, mat.cols);
    printf("A\n");
    printMatrix(mat);
    printf("\nB\n");
    printMatrix(vec1);
    doolitleLUP(mat, L, U, P);
    printf("\nL:\n");
    printMatrix(L);
    printf("\nU:\n");
    printMatrix(U);
    printf("\nP:\n");
    printMatrix(P);

    matrix vec2 = create_zero_matrix(vec1.rows, 1);
    solve_forward(L, vec1, vec2);
    printf("\nY:\n");
    printMatrix(vec2);
    solve_backward(U, vec2, vec1);
    printf("\nX'\n");
    printMatrix(vec1);
    transpose_inplace(P);
    matrix x = matmul(P, vec1);
    printf("\nX\n");
    printMatrix(x);
    printf("\nA*P^T*X\n");
    matrix AP = matmul(mat, P);
    matrix APX = matmul(AP, x);
    printMatrix(APX);
    destroy_matrix(APX);
    destroy_matrix(AP);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(mat);
    destroy_matrix(vec1);
    destroy_matrix(vec2);
    destroy_matrix(x);
}