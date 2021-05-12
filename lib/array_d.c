//
// Created by mateu on 11.05.2021.
//

#include "array_d.h"
#include <assert.h>
#include <stdio.h>

array create_array(size_t size) {
    double *data = malloc(size * sizeof(double));
    array arr = {size, data};
    return arr;
}

array create_zero_array(size_t size) {
    double *data = calloc(size, sizeof(double));
    array arr = {size, data};
    return arr;
}

void move_array(array *from, array *to) {
    free(to->data);
    *((size_t*) &to->size) = from->size;
    *((double **) &to->data) = from->data;
    *((size_t*) &from->size) = -1;
    *((double**) &from->data) = NULL;
}

void zero_array(array arr) {
    for (size_t i = 0; i < arr.size; ++i)
        arr.data[i] = 0.0;
}

double arrdot(array first, array second) {
    assert(first.size == second.size);
    double sum = 0.0;
    for (int i = 0; i < first.size; ++i) {
        sum += first.data[i] * second.data[i];
    }
    return sum;
}

void destroy_array(array arr) {
    if (arr.data != NULL)
        free(arr.data);
    *((size_t*) &arr.size) = -1;
    *((double**) &arr.data) = NULL;
}

void print_array(array arr) {
    for (int i = 0; i < arr.size; ++i) {
        printf("%lf ", arr.data[i]);
    }
    printf("\n");
}
