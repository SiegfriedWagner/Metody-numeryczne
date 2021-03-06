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
        print_matrix(result.output_matrix);
        matrix L = mat_create(result.output_matrix.rows, result.output_matrix.cols);
        matrix U = mat_create(result.output_matrix.rows, result.output_matrix.cols);
        matrix P = mat_create(result.output_matrix.rows, result.output_matrix.cols);
        zero_matrix(L);
        zero_matrix(U);
        zero_matrix(P);
        doolitleLUP(result.output_matrix, L, U, P);
        //doolitleLU(result.output_matrix, L, U);
        printf("L:\n");
        print_matrix(L);
        printf("U':\n");
        print_matrix(U);
        printf("P:\n");
        print_matrix(P);
        matrix recreated = mat_mul_mat(L, U);
        printf("L*U'");
        print_matrix(recreated);
        matrix LUP = mat_mul_mat(recreated, P);
        printf("L*U'*P\n");
        print_matrix(LUP);
        destroy_matrix(LUP);
        destroy_matrix(L);
        destroy_matrix(U);
    }
    fclose(file);
    destroy_matrix(result.output_matrix);
}