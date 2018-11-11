cmake_minimum_required (VERSION 3.12)
project (smlm LANGUAGES CXX CUDA)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(Matlab REQUIRED)


# -------------
# Mex files for CPU code
# -------------

matlab_add_mex(NAME mex_gausstransform_cpu SRC gausstransform/mex_gausstransform_cpu.cpp)

matlab_add_mex(NAME mex_expdist_cpu SRC expdist/mex_expdist_cpu.cpp)







option (USE_GPU
        "Use CUDA-enabled GPU Functions?" ON)
if (USE_GPU)


# -------------
# CUDA
# -------------

find_package(CUDA REQUIRED)
message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")

set(CMAKE_CUDA_FLAGS "-O3 -m64 -Xcompiler=-fPIC -Xptxas=-v -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60")
message(STATUS "CUDA compiler args cmake cuda flags ${CMAKE_CUDA_FLAGS}")



# -------------
# CUB
# -------------

find_package(CUB REQUIRED)
include_directories(${CUB_INCLUDE_DIR})

add_library(gausstransform SHARED gausstransform/gausstransform.cu)

add_library(expdist SHARED expdist/expdist.cu)


# Build GPU-enabled MEX functions



include_directories(${CUDA_INCLUDE_DIRS})
link_directories(${CUDA_LIBRARIES})

matlab_add_mex(NAME mex_gausstransform SRC gausstransform/mex_gausstransform.cpp LINK_TO gausstransform)

matlab_add_mex(NAME mex_expdist SRC expdist/mex_expdist.cpp LINK_TO expdist)





endif (USE_GPU)

