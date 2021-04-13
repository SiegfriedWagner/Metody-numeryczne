//
// Created by Mateusz Chojnowski on 13.04.2021.
//
#include "array.h"

array create_array(size_t size) {
    float *data = malloc(size * sizeof(float));
    array arr = {size, data};
    return arr;
}

array create_zero_array(size_t size) {
    float *data = calloc(size, sizeof(float));
    array arr = {size, data};
    return arr;
}

void move_array(array *from, array *to) {
    free(to->data);
    *((size_t*) &to->size) = from->size;
    *((float **) &to->data) = from->data;
    *((size_t*) &from->size) = -1;
    *((float**) &from->data) = NULL;
}

void destroy_array(array arr) {
    if (arr.data != NULL)
        free(arr.data);
    *((size_t*) &arr.size) = -1;
    *((float**) &arr.data) = NULL;
}
