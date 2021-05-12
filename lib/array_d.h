//
// Created by mateu on 11.05.2021.
//
#include <stdlib.h>

#ifndef NUMERYCZNE_ARRAY_D_H
#define NUMERYCZNE_ARRAY_D_H

typedef struct array {
    const size_t size;
    double * const data;
} array;

array create_array(size_t size);
array create_zero_array(size_t size);
void move_array(array *from, array *to);
void zero_array(array arr);
double arrdot(array first, array second);
void destroy_array(array arr);
void print_array(array arr);

#endif //NUMERYCZNE_ARRAY_D_H
