project(Zadanie_10)
add_executable(Zadanie_10 main.c)
add_executable(Oryginalne orginalny.cpp)
add_executable(Tests test.c)
target_link_libraries(Tests matrix_d array_d interop_d)
target_link_libraries(Zadanie_10 matrix_d array_d interop_d)