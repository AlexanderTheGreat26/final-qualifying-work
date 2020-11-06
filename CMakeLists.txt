cmake_minimum_required(VERSION 3.17)
project(ReflectionOnSquareBarriers)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")
add_executable(ReflectionOnSquareBarriers main.cpp)