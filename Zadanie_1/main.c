#include <stdio.h>
#include <float.h>
#include <string.h>
#include <libgen.h>

double sum_double_up(int n)
{
    double sum = 0;
    for (int i = 1; i <= n; ++i)
        sum += 1.0 / i;

    return sum;
}

double sum_double_down(int n)
{
    double sum = 0;
    for(int i = n; i >= 1; --i)
        sum += 1.0 / i;
    return sum;
}

float sum_float_up(int n)
{
    float sum = 0;
    for (int i = 1; i <= n; ++i)
        sum += 1.0f / i;
    return sum;
}

float sum_float_down(int n)
{
    float sum = 0;
    for (int i = n; i >= 1; --i)
        sum += 1.0f / i;
    return sum;
}

int main() {
    FILE *file;
    const int MAX_N = 1000000000;
    const int STEP_N = 10;
    int n = 1;
    float result_up, result_down;
    double dresult_up, dresult_down;
    char *path = strdup(__FILE__);
    path = dirname(path);
    strcat(path, "/sum.csv");
    file = fopen(path, "w");
    fprintf(file, "n,sum_float_up,sum_float_down,sum_double_up,sum_double_down\n");
    while(n < MAX_N)
    {
        result_up = sum_float_up(n);
        result_down = sum_float_down(n);
        dresult_up = sum_double_up(n);
        dresult_down = sum_double_down(n);
        fprintf(file, "%i,%.*f,%.*f,%.*f,%.*f\n", n, FLT_DIG, result_up, FLT_DIG, result_down, DBL_DIG, dresult_up, DBL_DIG, dresult_down);
        n *= STEP_N;
    }
    fclose(file);
    return 0;
}
