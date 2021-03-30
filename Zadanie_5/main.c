//
// Created by Mateusz Chojnowski on 24.03.2021.
//

#include "../lib/matrix.h"

int main(int argc, char *argv[]) {
    if(argc != 2) {
        printf("Missing argument: input file ");
        return 1;
    }
    unsigned int bufferSize = 40;
    char buffer[bufferSize];
    FILE *file = fopen(argv[1], "r");
    int eq_num = 0;
    if (readInt(file, &eq_num, buffer, bufferSize) != CORRECT){
        printf("Error while reading equations num");
        return 2;
    }
    matrix mat = create_matrix(eq_num, eq_num), L = create_zero_matrix(eq_num, eq_num), U = create_zero_matrix(eq_num, eq_num), P = create_zero_matrix(eq_num, eq_num);
    matrix bmat = create_matrix(eq_num, 1);
    // fill matrix
    for (int row = 0; row < mat.rows; ++row) {
        for (int col = 0; col < mat.cols; ++col) {
            float val = 0.0f;
            if (readFloat(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                return 2;
            }
            mat.data[row][col] = val;
        }
    }
    // fill b vector
    for (int row = 0; row < bmat.rows; ++row) {
        for (int col = 0; col < bmat.cols; ++col) {
            float val = 0.0f;
            if (readFloat(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                return 2;
            }
            bmat.data[row][col] = val;
        }
    }
    printMatrix(mat);
    printMatrix(bmat);
    doolitleLUP(mat, L, U, P);
    matrix y = create_zero_matrix(eq_num, 1);
    for (int row = 0; row < L.rows; ++row) {
        float b_val = bmat.data[row][0];
        for (int col = 0; col < row; ++col) {
            b_val -= L.data[row][col] * y.data[col][0];
        }
        y.data[row][0] = b_val / L.data[row][row];
    }
    matrix xprim = create_zero_matrix(eq_num, 1);
    for (int row = U.rows - 1; row >= 0 ; --row) {
        float y_val = y.data[row][0];
        for (int col = U.cols - 1; col > row; --col) {
            y_val -= U.data[row][col] * xprim.data[col][0];
        }
        xprim.data[row][0] = y_val / U.data[row][row];
    }
    printf("\n");
    transposeInplace(P);
    matrix x = matmul(P, xprim);
    printMatrix(x);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(mat);
    destroy_matrix(bmat);
}