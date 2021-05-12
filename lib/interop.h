//
// Created by mateu on 11.05.2021.
//

#ifndef NUMERYCZNE_INTEROP_H
#define NUMERYCZNE_INTEROP_H
#include "array.h"
#include "matrix.h"
void arr_to_mat_mul_h(matrix output, array first, array second);
matrix arr_to_mat_mul(array first, array second);
#endif //NUMERYCZNE_INTEROP_H
