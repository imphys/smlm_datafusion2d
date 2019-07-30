cmake_minimum_required (VERSION 3.14.3)
project (smlm LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR})
find_package(Matlab REQUIRED)


# -------------
# Mex files for CPU code
# -------------

matlab_add_mex(NAME mex_gausstransform_cpu SRC gausstransform/mex_gausstransform_cpu.cpp)

matlab_add_mex(NAME mex_expdist_cpu SRC expdist/mex_expdist_cpu.cpp)

install(TARGETS mex_gausstransform_cpu mex_expdist_cpu
        DESTINATION MATLAB/all2all)



# -------------
# Everything below this line is GPU related, using GPU is enabled by default
# To disable compilation of GPU code use cmake -DUSE_GPU=OFF .
# -------------


option (USE_GPU
        "Use CUDA-enabled GPU Functions?" ON)
if (USE_GPU)

#enable_language(CUDA)

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

message(STATUS "Found CUB ${CUB_INCLUDE_DIR}")

cuda_add_library(gausstransform SHARED gausstransform/gausstransform.cu)

cuda_add_library(expdist SHARED expdist/expdist.cu)


# Build GPU-enabled MEX functions



include_directories(${CUDA_INCLUDE_DIRS})
link_directories(${CUDA_LIBRARIES})

matlab_add_mex(NAME mex_gausstransform SRC gausstransform/mex_gausstransform.cpp LINK_TO gausstransform)

matlab_add_mex(NAME mex_expdist SRC expdist/mex_expdist.cpp LINK_TO expdist)

install(TARGETS mex_gausstransform mex_expdist
        DESTINATION MATLAB/all2all)




endif (USE_GPU)

