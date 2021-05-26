//
// Created by Mateusz Chojnowski on 18.05.2021.
//
#include "../lib/svd.h"
#include "../lib/matrix_d.h"
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
        code = readDouble(file, x.data+i, buffer, bufferSize);
        if (code != CORRECT) break;
        code = readDouble(file, y.data+i, buffer, bufferSize);
        if (code != CORRECT) break;
        code = readDouble(file, sigma.data+i, buffer, bufferSize);
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

double function(double x, double nth_function) {
    return pow(x, nth_function);
}

int main(int argc, char *argv[]) {
    if(argc != 3)
    {
        printf("Required arguments are: input file with observations and assumed degree of equation (polynomial");
        return 1;
    }
    int degree = strtol(argv[2], NULL, 10);
    array x = { 0, NULL }, y = { 0, NULL }, sigma = { 0, NULL };
    {
        FILE *file = fopen(argv[1], "r");
        parsing_code code = read_data_points(file, &x, &y, &sigma);
        fclose(file);
        if (code != CORRECT) {
            printf("Error while parsing file, error code %d", code);
            return 2;
        }
    }
    matrix A = mat_create(x.size, degree);
    matrix b = mat_create(y.size, 1);
    matrix a = mat_create_zero(A.rows, 1);
    for (int row = 0; row < A.rows; ++row) {
        for (int col = 0; col < A.cols; ++col)
            A.data[row][col] = function(x.data[row], (float) col) / sigma.data[row];
        b.data[row][0] = y.data[row] / sigma.data[row];
    }
    matrix Aorig = mat_create_copy(A);
    matrix U = mat_create_zero(A.rows, A.rows);
    matrix Sigma = mat_create_zero(A.rows, A.cols);
    matrix V = mat_create_zero(A.cols, A.cols);
    SVD(A, U, V, Sigma);
    array u_i = create_array(A.cols);
    array v_i = create_array(V.cols);
    // calculate a
    for (int i = 0; i < A.cols; ++i) {
        // copy u_i column
        for (int j = 0; j < A.cols; ++j)
            u_i.data[j] = U.data[j][i];
        printf("\nu_i(%d)\n");
        print_array(u_i);
        double U_ib = 0.0;
        // u_iT * b
        for (int j = 0; j < A.cols; ++j)
            U_ib += u_i.data[j] * b.data[j][0];
        for (int j = 0; j < v_i.size; ++j)
            v_i.data[j] = V.data[i][j] * U_ib / Sigma.data[i][i];
        // assert(a.rows == v_i.size);
        for (int j = 0; j < v_i.size; ++j)
            a.data[j][0] += v_i.data[j];
    }
    printf("\na(coefficients)\n");
    for (int i = 0; i < degree; ++i) {
        printf("%f ", a.data[i][0]);
    }
    printf("\n");
    destroy_matrix(b);
    destroy_matrix(A);
    destroy_array(x);
    destroy_array(y);
    destroy_array(sigma);
    return 0;
}
