cmake_minimum_required(VERSION 3.17)
project(numeryczne C)
set(CMAKE_C_STANDARD 99)

add_library(matrix STATIC lib/matrix.h lib/matrix.c)

add_subdirectory(Zadanie_1)
add_subdirectory(Zadanie_2)
add_subdirectory(Zadanie_3)
add_subdirectory(Zadanie_4)
add_subdirectory(Zadanie_5)