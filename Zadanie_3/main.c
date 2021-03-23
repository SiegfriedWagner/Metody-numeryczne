//
// Created by Mateusz Chojnowski on 10.03.2021.
//
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

void print_bytes(float element)
{
    unsigned int* as_int = &element;
    unsigned int filter = 1 << 31;
    int i = 0;
    printf("|");
    while (filter != 0)
    {
        if (i == 1 || i == 9)
            printf("|");
        i++;
        unsigned int var = ((*as_int) & filter);
        if (var)
            printf("%d", 1);
        else
            printf("%d", 0);
        filter /= 2;
    }
    printf("|\n");
}

void disect_float(float input, short* sign_out, short* exp_out, float* mantis_out) {
    unsigned int* as_int = &input;
    unsigned int filter = 1 << 31;
    if ((*as_int) & filter)
        *sign_out = -1;
    else
        *sign_out = 1;
    filter >>= 1;
    // exponent
    short reminder = 64;
    if (((*as_int) & filter))
        *exp_out = 1;
    else
        *exp_out = -127;
    filter >>= 1;
    for (int i = 0; i < 7; ++i) {
        if ((*as_int) & filter)
            *exp_out += reminder;
        reminder /= 2;
        filter >>= 1;
    }
    float mult = 0.5f;
    if(*exp_out == -127 && ((*as_int) & filter) == 0) { // denormalised
        (*exp_out)--;
        filter >>= 1;
        while((*as_int & filter) == 0 && filter != 0) {
            (*exp_out)--;
            filter /= 2;
        }
        mult = 1.0f;
        *mantis_out = 0.0f;
    }
    else {
        *mantis_out = 1.0f;
    }
    while(filter != 0) {
        if((*as_int) & filter)
            *mantis_out += mult;
        mult /= 2.0f;
        filter /= 2;
    }

}

int main(int argc, char* argv[]) {
    short sign = 0;
    short exp_out = 0;
    float mantis = 0;
    if (argc != 2) {
        printf("missing argument (float number)");
        return 1;
    }
    float number = strtod(argv[1], NULL);
    printf("Binary number representation: ");
    print_bytes(number);
    printf("Big endian\n"); // little is not supported
    disect_float(number, &sign, &exp_out, &mantis);
    if (sign > 0)
        printf("Sign: +\n");
    else
        printf("Sign: -\n");
    printf("Exponent (2^n): %d\n", exp_out);
    printf("Mantis (normalized): %f\n", mantis);
    printf("Recreated number: %.6e\n" , sign * powf(2, exp_out) * mantis);

    //TESTS
    disect_float(0.0f, &sign, &exp_out, &mantis);
    assert(0.0f == sign * powf(2, exp_out) * mantis);
    disect_float(3.75f, &sign, &exp_out, &mantis);
    assert(3.75f == powf(2, exp_out) * mantis);
    disect_float(4.0f, &sign, &exp_out, &mantis);
    assert(4.0f == powf(2, exp_out) * mantis);
    disect_float(12.0f, &sign, &exp_out, &mantis);
    assert(12.0f == powf(2, exp_out) * mantis);
    disect_float(-12.0f, &sign, &exp_out, &mantis);
    assert(-12.0f == sign * powf(2, exp_out) * mantis);
    disect_float(100000000000000000000000000000000000000.0f, &sign, &exp_out, &mantis);
    assert(100000000000000000000000000000000000000.0f == powf(2, exp_out) * mantis);
    disect_float(0.000005f, &sign, &exp_out, &mantis);
    assert(0.000005f == powf(2, exp_out) * mantis);
    disect_float(10e-41f, &sign, &exp_out, &mantis);
    assert(10e-41f == powf(2, exp_out) * mantis);
    return 0;
}