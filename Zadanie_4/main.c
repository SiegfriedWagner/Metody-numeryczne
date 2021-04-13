//
// Created by Mateusz Chojnowski on 17.03.2021.
//

#include <stdio.h>
#include "../lib/matrix.h"

int main(int argc, char **argv) {
    FILE* file;
    if(argc != 2) {
        printf("Missing argument: input file ");
        return 1;
    }

    file = fopen(argv[1], "r");
    parsing_result result = from_file(file);
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