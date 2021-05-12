//
// Created by Mateusz Chojnowski on 9.05.2021.
//
#include <math.h>
#include <assert.h>
#include <float.h>
#include "../lib/matrix_d.h"
#include "../lib/array_d.h"
#include "../lib/interop_d.h"

array diag(matrix mat) {
    size_t range = mat.rows > mat.cols ? mat.cols : mat.rows;
    array returned = create_array(range);
    for (int i = 0; i < range; ++i)
        returned.data[i] = mat.data[i][i];
    return returned;
}

array not_diag(matrix mat) {
    size_t range = mat.rows > mat.cols ? mat.cols : mat.rows;
    array returned = create_array(range);
    for (int i = 1; i < range; ++i)
        returned.data[i] = mat.data[i][i - 1];
    returned.data[0] = 0.0;
    return returned;
}

int sign(double num) {
    return (0 < num) - (num < 0);
}

void calculate_vector_x_k(array x, array alfa, matrix A_k, size_t k) {
    // calculate s_k
    double s_k = 0.0;
    for (int i = 0; i < A_k.rows; ++i)
        s_k += pow(A_k.data[i][k], 2.0);
    s_k = sqrt(s_k);
    // a_k = -s_k * sign(A_kk)
    double alfa_k = -s_k * sign(A_k.data[k][k]);
    alfa.data[k] = alfa_k;
    // vector x
    if (s_k == 0) {
        zero_array(x);
        alfa.data[k] = 0.0f;
    } else {
        double c_k = 0.0f;
        for (int i = 0; i < A_k.rows; ++i) {
            if (i < k)
                x.data[i] = 0.0f;
            else if (i == k) {
                x.data[i] = sqrt(0.5 * (1 + (fabs(A_k.data[k][k]) / s_k)));
                c_k = sign(A_k.data[k][k]) / (2.0f * s_k * x.data[k]);
            } else
                x.data[i] = c_k * A_k.data[i][k];
        }
    }
}

void calculate_vector_y_k(array y, array beta, matrix A_k, size_t k) {
    double t_k = 0.0f;
    for (size_t j = k + 1; j < A_k.cols; ++j) {
        t_k += pow(fabs(A_k.data[k][j]), 2.0);
    }
    t_k = sqrt(t_k);
    double b_k = -t_k * sign(A_k.data[k][k + 1]);
    beta.data[k] = b_k;
    if (t_k == 0.0) {
        zero_array(y);
        beta.data[k] = 0.0;
    } else {
        double d_k = 0.0;
        for (int j = 0; j < A_k.cols; ++j) {
            if (j <= k)
                y.data[j] = 0.0;
            else if (j == k + 1) {
                y.data[j] = sqrt(0.5 * (1 + fabs(A_k.data[k][j]) / t_k));
            } else {
                d_k = ((double) sign(A_k.data[k][j])) / (2.0 * t_k * y.data[j]);
                y.data[j] = d_k * A_k.data[k][j];
            }
        }
    }
}

void calculate_P_or_Q(matrix P_or_Q, array x_or_y) {
    assert(P_or_Q.rows == P_or_Q.cols);
    matrix unit = mat_create_identity(P_or_Q.rows);
    matrix X = arr_to_mat_mul(x_or_y, x_or_y);
    mat_mul_scalar_inplace(X, 2.0);
    mat_sub_mat_h(unit, X, P_or_Q);
    destroy_matrix(X);
    destroy_matrix(unit);
//    arr_to_mat_mul_h(P_or_Q, x_or_y, x_or_y);
//    mat_mul_scalar_inplace(P_or_Q, -2.0);
//    for (int i = 0; i < P_or_Q.rows; ++i)
//        P_or_Q.data[i][i] += 1.0;
}

int sturm(double upper, array diagValues, array notDiag) {
    array q = create_zero_array(diagValues.size);
    int returned = 0;
    q.data[0] = diagValues.data[0] - upper;
    for (int i = 1; i < q.size; i++) {

        if (fabs(q.data[i - 1]) > DBL_EPSILON)
            q.data[i] = (diagValues.data[i] - upper) - (notDiag.data[i] * notDiag.data[i]) / q.data[i - 1];
        else {
            q.data[i] = (diagValues.data[i] - upper) - fabs(notDiag.data[i]) / DBL_EPSILON;
        }
        if (q.data[i] < 0)
            returned++;
    }
    destroy_array(q);
    return returned;
}

array bisection(array diagValues, array notDiag, double lambda_min, double lambda_max) {
    array top = create_array(diagValues.size);
    array bottom = create_array(diagValues.size);
    for (int i = 0; i < top.size; ++i) {
        top.data[i] = lambda_max;
        bottom.data[i] = lambda_min;
    }
    for (int k = 0; k < top.size - 1; --k) {
        double upper = top.data[k];
        double lower = bottom.data[k];
        double value = 0.0;
        while (fabs(upper - lower) > 0) {
            if (value == upper - lower)
                break;
            else {
                value = upper - lower;
                top.data[k] = (upper - lower) / 2.0;
                int a = sturm(top.data[k], diagValues, notDiag);
                if (a <= k)
                    lower = top.data[k];
                else
                    upper = top.data[k];
                for (int i = 0; i < a; i++)
                    if (top.data[i] > top.data[k])
                        top.data[i] = top.data[k];
                for (int i = a; i < k; i++)
                    if (bottom.data[i] < top.data[k])
                        bottom.data[i] = top.data[k];
            }
        }
    }
    destroy_array(bottom);
    return top;
}

void calculate_J(matrix A, matrix P, matrix Q, matrix J) {
    assert(A.rows >= A.cols);
    assert(P.rows == A.rows && J.rows == A.rows && Q.rows == A.cols);
    assert(P.cols == A.rows && J.cols == A.cols && Q.cols == A.cols);
    int rows = A.rows;
    int cols = A.cols;
    matrix tempA = mat_create_zero(rows, cols);
    matrix Qk = mat_create_zero(cols, cols);
    matrix Qkm1 = mat_create_identity(cols);
    matrix Pk = mat_create_zero(rows, rows);
    matrix Pkm1 = mat_create_identity(rows);
    array x = create_zero_array(rows);
    array y = create_zero_array(cols);
    array alfa = create_zero_array(rows);
    array beta = create_zero_array(cols);

    for (int k = 0; k < cols; ++k) {
        // calculate x_k
        calculate_vector_x_k(x, alfa, A, k);
        //printf("\nalfa\n");
        //print_array(alfa);
        // calculate P
        calculate_P_or_Q(P, x);
        mat_mul_mat_h(Pkm1, P, Pk);
        copy_values(Pk, Pkm1);
        mat_mul_mat_h(P, A, tempA);
//        printf("\nP\n");
//        print_matrix(P);
//        printf("\nA\n");
//        print_matrix(A);
//        printf("\ntempA\n");
//        print_matrix(tempA);
        zero_matrix(P);
        if (k == cols - 1) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++)

                    J.data[i][j] = tempA.data[i][j];
            }
        }
        calculate_vector_y_k(y, beta, tempA, k);
        //printf("\ny\n");
        //print_array(y);
        //printf("\nbeta\n");
        //print_array(beta);
        // calculate Q
        calculate_P_or_Q(Q, y);
        mat_mul_mat_h(Qkm1, Q, Qk);
        copy_values(Qk, Qkm1);
        mat_mul_mat_h(tempA, Q, A);
        zero_matrix(Q);
    }
    copy_values(Qk, Q);
    copy_values(Pk, P);
    destroy_array(beta);
    destroy_array(alfa);
    destroy_array(y);
    destroy_array(x);
    destroy_matrix(Pkm1);
    destroy_matrix(Pk);
    destroy_matrix(Qkm1);
    destroy_matrix(Qk);
    destroy_matrix(tempA);
}

matrix calculate_Y(matrix K, array eigenK) {
    assert(K.rows == K.cols);
    matrix temp = mat_create(K.rows, K.cols);
    matrix Y = mat_create_zero(K.rows, K.cols);
    matrix L = mat_create(K.rows, K.cols), U = mat_create(K.rows, K.cols), P = mat_create(K.rows, K.cols);
    matrix A = mat_create_zero(K.rows, K.cols);
    matrix solution = mat_create_zero(K.rows, 1);
    matrix col = mat_create_zero(K.rows, 1);
    for (int i = 0; i < A.cols; ++i) {
        diagonal(temp, 1.0);
        mat_mul_scalar_inplace(temp, eigenK.data[eigenK.size - 1 - i]);
        mat_sub_mat_h(K, temp, A);
        doolitleLUP(A, L, U, P);
        solution.data[solution.rows - 1][0] = 1.0f;
        solve_backward_offset(U, solution, solution, 1);
        print_matrix(solution);
        mat_mul_mat_h(P, solution, col);
        mat_mul_scalar_inplace(col, 1.0 / mat_sum(col)); // TODO: Add divde function
        copy_column(col, Y, 0, i);
    }
    destroy_matrix(L);
    destroy_matrix(U);
    destroy_matrix(P);
    destroy_matrix(A);
    destroy_matrix(solution);
    destroy_matrix(col);
    return Y;
}

matrix calculate_X(matrix Y, matrix J, array szcz) {
    assert(Y.rows == Y.cols);
    size_t n = J.cols;
    size_t m = J.rows;
    matrix vec1 = mat_create_zero(n, 1);
    matrix vec2 = mat_create_zero(m, 1);
    matrix X = mat_create_zero(J.rows, J.cols);
    for (int col = 0; col < Y.cols; ++col) {
        copy_column(Y, vec1, col, 0);

        mat_mul_mat_h(J, vec1, vec2);
        double value = 1 / szcz.data[col];
        for (int i = 0; i < n; ++i) {
            X.data[i][col] = value * vec2.data[i][0];
        }
    }
    destroy_matrix(vec1);
    destroy_matrix(vec2);
    return X;
}

void SVD(matrix A, matrix U, matrix V, matrix Sigma) {
    assert(U.rows == A.rows && Sigma.rows == A.rows && V.rows == A.cols);
    assert(U.cols == A.rows && Sigma.cols == A.cols && V.cols == A.cols);
    matrix P = mat_create(A.rows, A.rows);
    matrix Q = mat_create(A.cols, A.cols);
    matrix J = mat_create(A.rows, A.cols);
    calculate_J(A, P, Q, J);
    printf("\nP\n");
    print_matrix(P);
    printf("\nJ\n");
    print_matrix(J);
    printf("\nQ\n");
    print_matrix(Q);
    matrix JT = transpose(J);
    matrix QT = transpose(Q);
    matrix PJ = mat_mul_mat(P, J);
    matrix PJQ = mat_mul_mat(PJ, QT);
    printf("\nPJQ\n");
    print_matrix(PJQ);
    copy_values(PJQ, A);
    matrix K = mat_mul_mat(JT, J);
    printf("\nK\n");
    // diag
    array diag_values = diag(K);
    array not_diag_values = not_diag(K);
    // lambda
    double lambda_min = 0.0;
    double lambda_max = 0.0;
    for (int i = 0; i < diag_values.size - 1; i++) {
        double min_value = diag_values.data[i] - (fabs(not_diag_values.data[i]) + fabs(not_diag_values.data[i + 1]));
        if (lambda_min > min_value)
            lambda_min = min_value;
        double max_value = diag_values.data[i] + (fabs(not_diag_values.data[i]) + fabs(not_diag_values.data[i + 1]));
        if (lambda_max < max_value)
            lambda_max = max_value;
    }
    printf("\ndiag\n");
    print_array(diag_values);
    printf("\nnot diag\n");
    print_array(not_diag_values);
    array eigenK = bisection(diag_values, not_diag_values, lambda_min, lambda_max);
    // Y matrix
    matrix Y = calculate_Y(K, eigenK);
    printf("\nY\n");
    print_matrix(Y);
    // X matrix
    array szcz = create_array(K.cols);
    for (int i = 0; i < K.cols; ++i)
        szcz.data[i] = sqrt(eigenK.data[eigenK.size - 1 - i]);
    matrix X = calculate_X(Y, J, szcz);
    printf("\nX\n");
    print_matrix(Y);
    // Sigma matrix
    zero_matrix(Sigma);
    for (int i = 0; i < Sigma.rows; ++i) {
        Sigma.data[i][i] = szcz.data[i];
    }
    printf("\nSigma\n");
    print_matrix(Sigma);
    // U and V
    mat_mul_mat_h(P, X, U);
    printf("\nU\n");
    print_matrix(U);
    mat_mul_mat_h(Q, Y, V);
    printf("\nV\n");
    print_matrix(V);
    matrix VT = transpose(V);
    matrix temp = mat_mul_mat(U, Sigma);
    matrix temp2 = mat_mul_mat(temp, VT);
    printf("\nU*Sigma*V\n");
    print_matrix(temp2);
    destroy_matrix(temp2);
    destroy_matrix(temp);
    destroy_matrix(VT);
    destroy_matrix(X);
    destroy_matrix(Y);
    destroy_array(not_diag_values);
    destroy_array(diag_values);
    destroy_matrix(K);
    destroy_matrix(PJ);
    destroy_matrix(QT);
    destroy_matrix(JT);
    destroy_matrix(J);
    destroy_matrix(Q);
    destroy_matrix(P);
}

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
    SVD(A, U, V, Sigma);
    destroy_matrix(U);
    destroy_matrix(Sigma);
    destroy_matrix(V);
    destroy_matrix(A);
    return 0;
}