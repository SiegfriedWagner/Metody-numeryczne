//
// Created by Mateusz Chojnowski on 07.04.2021.
//
#include "../lib/matrix.h"
#include "../lib/array.h"
#include <math.h>
extern unsigned short bufferSize;


parsing_code read_data_points(FILE *file, array *x_out, array *y_out, array *sigma_out) {
    char buffer[bufferSize];
    int size = 0;
    parsing_code code = readInt(file, &size, buffer, bufferSize);
    if (code != CORRECT)
        return code;
    array x = create_array(size);
    array y = create_array(size);
    array sigma = create_array(size);
    for (int i = 0; i < size; ++i) {
        code = readFloat(file, x.data+i, buffer, bufferSize);
        if (code != CORRECT) break;
        code = readFloat(file, y.data+i, buffer, bufferSize);
        if (code != CORRECT) break;
        code = readFloat(file, sigma.data+i, buffer, bufferSize);
        if (code != CORRECT) break;
    }
    if(code != CORRECT)
    {
        destroy_array(x);
        destroy_array(y);
        destroy_array(sigma);
    }
    else {
        move_array(&x, x_out);
        move_array(&y, y_out);
        move_array(&sigma, sigma_out);
    }
    return code;
}

float function(float x, float nth_function) {
    return powf(x, nth_function);
}

int main(int argc, char* argv[]) {

    if(argc != 3)
    {
        printf("Required arguments are: input file with observations and assumed degree of equation (polynomial");
        return 1;
    }
    int degree = strtol(argv[2], NULL, 10);
    array x = { -1, NULL }, y = { -1, NULL}, sigma = { -1, NULL};
    {
        FILE *file = fopen(argv[1], "r");
        parsing_code code = read_data_points(file, &x, &y, &sigma);
        fclose(file);
        if (code != CORRECT) {
            printf("Error while parsing file, error code %d", code);
            return 2;
        }
    }
    matrix A = create_matrix(x.size, degree);
    matrix b = create_matrix(y.size, 1);
    for (int row = 0; row < A.rows; ++row) {
        for (int col = 0; col < A.cols; ++col) {
            A.data[row][col] = function(x.data[row], (float) col) / sigma.data[row];
        }
        b.data[row][0] = (y.data[row] ) / sigma.data[row];
    }
    printf("\nA\n");
    printMatrix(A);
    printf("\nb\n");
    printMatrix(b);
    matrix AT = transpose(A);
    matrix alpha = matmul(AT, A);
    printf("\nalpha\n");
    printMatrix(alpha);
    matrix alpha_inv = copy_matrix(alpha);
    matrix beta = matmul(AT, b);
    matrix a = solve_equation(alpha, beta);
    printf("\ncoefficients (from lowest to highest)\n");
    printMatrix(a);
    invert_matrix_inplace(alpha_inv);
    printf("\ncond(alpha) = %f\n", norm(alpha) * norm(alpha_inv));
    float chi = 0;
    // calculate chi^2
    for (int i = 0; i < y.size; ++i) {
        float fun_sum = 0.0f;
        for (int d = 0; d < degree; ++d) {
            fun_sum += a.data[d][0] * function(x.data[i], (float) d);
        }
        float temp = (y.data[i] - fun_sum) / sigma.data[i];
        chi += temp * temp;
    }
    printf("Chi^2 = %f", chi);
    printf("\nCov matrix - alpha^(-1)\n");
    printMatrix(alpha_inv);
    destroy_array(x);
    destroy_array(y);
    destroy_array(sigma);
    destroy_matrix(alpha_inv);
    destroy_matrix(a);
    destroy_matrix(beta);
    destroy_matrix(alpha);
    destroy_matrix(AT);
    destroy_matrix(b);
    destroy_matrix(A);
}
