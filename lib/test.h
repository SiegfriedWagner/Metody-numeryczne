//
// Created by mateu on 18.05.2021.
//

#ifndef NUMERYCZNE_TEST_H
#define NUMERYCZNE_TEST_H

#include <stdio.h>
#include "math.h"
#define absolute(N) (((N)<0)?(-(N)):(N)) \

#define all_equal1(truth, tested, size, delta_max)                                                                     \
({                                                                                                                     \
   int equal = 1;                                                                                                      \
   for(int i = 0; i < (size); i++) {                                                                                   \
      if (absolute((truth)[i] - (tested)[i]) > (delta_max)) {                                                          \
         equal = 0;                                                                                                    \
         printf("values at index[%d] are unequal\ntruth: %f\n tested: %f\n", i,  (truth)[i], (tested)[i]);             \
      }                                                                                                                \
   }                                                                                                                   \
   equal;                                                                                                              \
})


#define all_equal2(truth, tested, rows, columns, delta_max)                                                            \
({                                                                                                                     \
   bool equal = true;                                                                                                  \
   for(int i = 0; i < (rows); i++)                                                                                     \
      for(int j = 0; j < (columns); j++)                                                                               \
         if (absolute((truth)[i][j] - (tested)[i][j]) > (delta_max)) {                                                      \
            equal = false;                                                                                             \
            printf("values at index[%d, %d] are unequal\ntruth: %f\n tested: %f\n", i, j, (truth)[i][j], (tested)[i][j]); \
         }                                                                                                             \
   equal;                                                                                                              \
})

#endif //NUMERYCZNE_TEST_H
