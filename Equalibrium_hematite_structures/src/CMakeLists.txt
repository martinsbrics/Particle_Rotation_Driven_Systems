cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_FLAGS "-fopenmp -ffast-math  -Wall  -Wno-write-strings -Wno-unused-variable -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-value -Wno-unused-result -std=c++11")#
#set(CMAKE_CXX_FLAGS "  -Wall  -Wno-write-strings -Wno-unused-variable -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-value -Wno-unused-result -std=c++11")#
LINK_DIRECTORIES( ${LINK_DIRECTORIES} ~/bin/lib /usr/local/lib) 
include_directories(${include_directories} ~/bin/include/eigen3 /usr/local/include/optim) 

project(two_cube_2d_new)

add_executable(two_cube_2d_new main.cpp calculate.cpp calc.cpp)
 
TARGET_LINK_LIBRARIES(two_cube_2d_new lapack blas fftw3 fftw3_omp m  rt optim gsl)#gfortran umfpack amd superlu 

install(TARGETS two_cube_2d_new RUNTIME DESTINATION bin)
