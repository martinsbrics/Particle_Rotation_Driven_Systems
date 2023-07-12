cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_FLAGS "-fopenmp -ffast-math  -Wall  -Wno-write-strings -Wno-unused-variable -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-value -Wno-unused-result -std=c++11")#
#set(CMAKE_CXX_FLAGS "  -Wall  -Wno-write-strings -Wno-unused-variable -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-value -Wno-unused-result -std=c++11")#
LINK_DIRECTORIES( ${LINK_DIRECTORIES} ~/bin/lib /usr/local/lib) 
include_directories(${include_directories} ~/bin/include/eigen3 /usr/local/include/optim) 

project(micrascopic_filament)

add_executable(micrascopic_filament main.cpp integrator.cpp simulation.cpp)

install(TARGETS micrascopic_filament RUNTIME DESTINATION bin)
