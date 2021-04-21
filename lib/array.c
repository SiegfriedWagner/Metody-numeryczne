//
// Created by Mateusz Chojnowski on 13.04.2021.
//
#include "array.h"
#include <assert.h>
#include <stdio.h>

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

float arrdot(array first, array second) {
    assert(first.size == second.size);
    float sum = 0.0f;
    for (int i = 0; i < first.size; ++i) {
        sum += first.data[i] * second.data[i];
    }
    return sum;
}

void destroy_array(array arr) {
    if (arr.data != NULL)
        free(arr.data);
    *((size_t*) &arr.size) = -1;
    *((float**) &arr.data) = NULL;
}

void print_array(array arr) {
    for (int i = 0; i < arr.size; ++i) {
        printf("%f ", arr.data[i]);
    }
    printf("\n");
}
