//
// Created by mateu on 03.03.2021.
//

#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#define FLT_EPSILON_FULL_PRINT 23
#define FLT_MIN_FULL_PRINT 54
#define FLT_LOWEST_FULL_PRINT 61
#define DBL_EPSILON_FULL_PRINT 32


float calculate_emach(FILE* log_file, int log_stdout, int* mantysa)
{
    // calculates e_mach and logs each step to file
    if (log_stdout)
        printf("n,emach,1+emach\n");
    if (log_file != NULL)
        fprintf(log_file, "n, emach,1+emach\n");

    *mantysa = 1; // zero bit
    float emach = 1.0f;
    float newemach = emach / 2.0f;
    while(1.0f + newemach > 1.0f)
    {
        emach = newemach;
        if(log_stdout)
            printf("%d,%.*f,%.*f\n", *mantysa, FLT_EPSILON_FULL_PRINT, emach, FLT_EPSILON_FULL_PRINT, 1.0f + emach);
        if(log_file != NULL)
            fprintf(log_file, "%d,%.*f,%.*f\n", *mantysa, FLT_EPSILON_FULL_PRINT, emach, FLT_EPSILON_FULL_PRINT, 1.0f + emach);
        newemach /= 2.0f;
        (*mantysa)++;
    }
    if(log_stdout)
        printf("%d,%.*f,%.*f\n", ++(*mantysa), FLT_EPSILON_FULL_PRINT, newemach, FLT_EPSILON_FULL_PRINT, 1.0f + newemach);
    if(log_file != NULL)
        fprintf(log_file, "%d,%.*f,%.*f\n", (*mantysa), FLT_EPSILON_FULL_PRINT, newemach, FLT_EPSILON_FULL_PRINT, 1.0f + newemach);
    return emach;
}


float calculate_lowest_positive_float()
{
    float ret = 1.0f;
    float num = ret / 2.0f;
    while (num)
    {
        //printf("test\n");
        ret = num;
        num /= 2.0f;
    }
    return ret;
}

float calculate_min_positive_float()
{
    int _;
    // calcuates minumum positive float in normalized form
    float lowest = 1.0f +  calculate_emach(NULL, 0, &_);
    float newlowest = lowest / 2.0f;
    while (lowest == (newlowest * 2.0f))
    {
        lowest = newlowest;
        newlowest /= 2.0f;

    }
    return newlowest * 2.0f;
}

float calculate_highest_float()
{
    float highest = 1.0f;
    float newemach = 1.0f / 2.0f;
    while(1.0f + newemach > 1.0f)
    {
        highest += newemach;
        newemach = newemach / 2.0f;
    }
    float newhighest = highest * 2.0f;
    while(highest == newhighest / 2.0f)
    {
        highest = newhighest;
        newhighest *= 2.0f;
    }
    return highest;
}

double calculate_emach_double(FILE* log_file, int log_stdout, int* mantysa)
{
    // calculates double e_mach and logs each step to file
    if (log_stdout)
        printf("n,emach,1+emach\n");
    if (log_file != NULL)
        fprintf(log_file, "n,emach,1+emach\n");
    *mantysa = 1; // zero bit
    double emach = 1.0f;
    double newemach = emach / 2.0;
    while(1.0f + newemach > 1.0f)
    {

        emach = newemach;
        if(log_stdout)
            printf("%d,%.*f,%.*f\n", (*mantysa), DBL_EPSILON_FULL_PRINT, emach, DBL_EPSILON_FULL_PRINT, 1.0 + emach);
        if(log_file != NULL)
            fprintf(log_file, "%d,%.*f,%.*f\n", (*mantysa), DBL_EPSILON_FULL_PRINT, emach, DBL_EPSILON_FULL_PRINT, 1.0 + emach);
        emach = newemach;
        newemach /= 2.0;
        (*mantysa)++;
    }
    if(log_stdout)
        printf("%d,%.*f,%.*f\n", (*mantysa), DBL_EPSILON_FULL_PRINT, newemach, DBL_EPSILON_FULL_PRINT, 1.0 + newemach);
    if(log_file != NULL)
        fprintf(log_file, "%d,%.*f,%.*f\n", (*mantysa), DBL_EPSILON_FULL_PRINT, newemach, DBL_EPSILON_FULL_PRINT, 1.0 + newemach);
    return emach;
}

double calculate_lowest_positive_double()
{
    double ret = 1.0;
    double num = ret / 2.0;
    while (num)
    {
        //printf("test\n");
        ret = num;
        num /= 2.0;
    }
    return ret;
}

double calculate_min_positive_double()
{
    int _;
    // calcuates minumum positive double in normalized form
    double lowest = 1.0 + calculate_emach_double(NULL, 0, &_);
    double newlowest = lowest / 2.0;
    while (lowest == (newlowest * 2.0))
    {
        lowest = newlowest;
        newlowest /= 2.0;
    }
    return newlowest * 2.0;
}

double calculate_highest_double()
{
    double highest = 1.0;
    double newemach = 1.0 / 2.0;
    while(1.0 + newemach > 1.0)
    {
        highest += newemach;
        newemach = newemach / 2.0;
    }
    double newhighest = highest * 2.0f;
    while(highest == newhighest / 2.0f)
    {
        highest = newhighest;
        newhighest *= 2.0f;
    }
    return highest;
}


int main() {
    FILE *file;
    file = fopen("emach_float.csv", "w");
    int mantis = 0;
    float float_emach = calculate_emach(file, 0, &mantis);
    fclose(file);
    printf("e_mach:                    %.6e (%.*f)\n", float_emach, FLT_EPSILON_FULL_PRINT, float_emach);
    assert(FLT_EPSILON == float_emach);
    printf("float_mantis:              %d\n", mantis);
    assert(mantis == FLT_MANT_DIG);
    float float_min_normalized = calculate_min_positive_float();
    assert(float_min_normalized == FLT_MIN);
    printf("float_min_normalized:      %.6e (%.*f)\n", float_min_normalized, FLT_MIN_FULL_PRINT, float_min_normalized);
    float float_min_denormalized = calculate_lowest_positive_float();
    printf("float_min_denormalized:    %.6e (%.*f)\n",float_min_denormalized, FLT_LOWEST_FULL_PRINT, float_min_denormalized);
    assert(float_min_denormalized == pow(2, -149));
    float float_max = calculate_highest_float();
    printf("float_max:                 %.6e (%.f)\n", float_max, float_max);
    assert(float_max == FLT_MAX);
    file = fopen("emach_double.csv", "w");
    double double_emach = calculate_emach_double(file, 0, &mantis);
    fclose(file);
    printf("e_mach_double:             %.*f\n", DBL_EPSILON_FULL_PRINT, double_emach);
    assert(double_emach == DBL_EPSILON);
    printf("double_mantis              %.d\n", mantis);
    assert(mantis == DBL_MANT_DIG);
    double double_min_normalized = calculate_min_positive_double();
    printf("double_min_normalized:     %.6e\n", double_min_normalized);
    assert(double_min_normalized == DBL_MIN);
    double double_min_denormalized = calculate_lowest_positive_double();
    printf("double_min_denormalized:   %.6e\n", double_min_denormalized);
    assert(double_min_denormalized == pow(2, -(1022 + 52)));
    double double_max = calculate_highest_double();
    printf("double_max:                %.6e\n", double_max);
    assert(double_max == DBL_MAX);
    return 0;

}