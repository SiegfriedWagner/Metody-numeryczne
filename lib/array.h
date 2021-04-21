//
// Created by mateu on 13.04.2021.
//
#include <stdlib.h>

#ifndef ARRAY_H
#define ARRAY_H

typedef struct array {
    const size_t size;
    float * const data;
} array;

array create_array(size_t size);
array create_zero_array(size_t size);
void move_array(array *from, array *to);
float arrdot(array first, array second);
void destroy_array(array arr);
void print_array(array arr);
#endif //NUMERYCZNE_VECTOR_H
