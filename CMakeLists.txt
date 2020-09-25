cmake_minimum_required (VERSION 3.9)

#>> Add in CUDA toolkit

#>> Set up Project and Subdirectories
project (strainlines LANGUAGES CXX)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set (CMAKE_CXX_COMPILER icpc)

#>> Find the C++ required libraries

find_package(OpenMP REQUIRED)

#>> Add in the source files

add_executable(strainlines
    src/main.cpp
)

#>> Set CUDA flags

#target_compile_options(strainlines PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O2 -arch sm_61 -Xcompiler "-fopenmp -mtune=native --fast-math" >)
target_compile_options(strainlines PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-qopenmp -ax=AVX -Wall -O3 -mkl>)# -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS -DARMA_NO_DEBUG -mkl>)

#>> Add includes

target_include_directories(strainlines PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/usr/local/include>)

#>> Add linkers

target_link_libraries(strainlines PRIVATE "${OpenMP_CXX_FLAGS} -mkl")
