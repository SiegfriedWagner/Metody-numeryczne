//
// Created by Mateusz Chojnowski on 11.05.2021.
//

#include "interop.h"
#include <assert.h>

void arr_to_mat_mul_h(matrix output, array first, array second) {
    assert(output.rows == first.size);
    assert(output.cols == second.size);
    for (int i = 0; i < first.size; ++i) {
        for (int j = 0; j < second.size; ++j) {
            output.data[i][j] = first.data[i] * second.data[j];
        }
    }
}

matrix arr_to_mat_mul(array first, array second) {
    matrix result = mat_create(first.size, second.size);
    arr_to_mat_mul_h(result, first, second);
    return result;
}
