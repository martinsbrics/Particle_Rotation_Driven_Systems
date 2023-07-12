cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS " -ffast-math  -Wall  -Wno-write-strings -Wno-unused-variable -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-value -Wno-unused-result -std=c++11")#-fopenmp
project(filament)
include_directories(${include_directories} ~/bin/include/eigen3) 

add_executable(filament main.cpp integrator.cpp integrator1.cpp simulation.cpp)

TARGET_LINK_LIBRARIES(filament   blas lapack)

install(TARGETS filament RUNTIME DESTINATION bin)
