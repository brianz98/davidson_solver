cmake_minimum_required (VERSION 3.12)
project (davidson)
enable_language (Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set (dialect "-ffree-form -std=f2008 -fimplicit-none -fopenmp -O3 -ffast-math")
    set (extras "-fcheck=all -Wall -ggdb -O0")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set (dialect "-stand f08 -free -implicitnone -qopenmp -O3")
    set (extras "-check bounds -O0")
endif ()

set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${extras}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

file (GLOB_RECURSE sources src/*.f90 src/*.F90)
add_executable(davidson.x ${sources})
target_link_libraries (davidson.x /opt/OpenBLAS/lib/libopenblas.a -lpthread)