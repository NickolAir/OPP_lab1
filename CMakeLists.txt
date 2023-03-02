cmake_minimum_required(VERSION 3.24)
project(SLAE C)

set(CMAKE_C_STANDARD 99)
include_directories(/usr/include/x86_64-linux-gnu/mpich)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

add_executable(SLAE main.c ../SLAE/main.c)
target_link_libraries(SLAE m)