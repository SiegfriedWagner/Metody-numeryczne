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
    print_matrix(mat);
    printf("\nB\n");
    print_matrix(vec1);
    doolitleLUP(mat, L, U, P);
    printf("\nL:\n");
    print_matrix(L);
    printf("\nU:\n");
    print_matrix(U);
    printf("\nP:\n");
    print_matrix(P);

    matrix vec2 = create_zero_matrix(vec1.rows, 1);
    solve_forward(L, vec1, vec2);
    printf("\nY:\n");
    print_matrix(vec2);
    solve_backward(U, vec2, vec1);
    printf("\nX'\n");
    print_matrix(vec1);
    transpose_inplace(P);
    matrix x = matmul(P, vec1);
    printf("\nX\n");
    print_matrix(x);
    printf("\nA*P^T*X\n");
    matrix AP = matmul(mat, P);
    matrix APX = matmul(AP, x);
    print_matrix(APX);
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