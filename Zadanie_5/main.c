//
// Created by Mateusz Chojnowski on 24.03.2021.
//

#include "../lib/matrix.h"
void solve_forward(matrix left, matrix right, matrix output) {
    assert(right.cols == output.cols && right.rows == output.rows);
    for (int row = 0; row < left.rows; ++row) {
        float b_val = right.data[row][0];
        for (int col = 0; col < row; ++col) {
            b_val -= left.data[row][col] * output.data[col][0];
        }
        output.data[row][0] = b_val / left.data[row][row];
    }
}

void solve_backward(matrix left, matrix right, matrix output) {
    assert(right.cols == output.cols && right.rows == output.rows);
    for (int row = left.rows - 1; row >= 0 ; --row) {
        float y_val = right.data[row][0];
        for (int col = left.cols - 1; col > row; --col) {
            y_val -= left.data[row][col] * output.data[col][0];
        }
        output.data[row][0] = y_val / left.data[row][row];
    }
}

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
    matrix vec1 = create_matrix(eq_num, 1);
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
    for (int row = 0; row < vec1.rows; ++row) {
        for (int col = 0; col < vec1.cols; ++col) {
            float val = 0.0f;
            if (readFloat(file, &val, buffer, bufferSize) != CORRECT) {
                printf("Error during parsing matrix at (%i, %i)", row, col);
                return 2;
            }
            vec1.data[row][col] = val;
        }
    }
    printMatrix(mat);
    printMatrix(vec1);
    doolitleLUP(mat, L, U, P);
    matrix vec2 = create_zero_matrix(eq_num, 1);
    solve_forward(L, vec1, vec2);
    solve_backward(U, vec2, vec1);
    printf("\n");
    transposeInplace(P);
    matrix x = matmul(P, vec1);
    printMatrix(x);
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(mat);
    destroy_matrix(vec1);
    destroy_matrix(vec2);
    destroy_matrix(x);
}