#!/usr/bin/env python

import numpy

from kernel_tuner import tune_kernel

from tune_utils import get_kernel_path


def tune_gausstransform():

    #setup test input
    size = numpy.int32(5000)
    ndim = numpy.int32(2)
    A = numpy.random.randn(size*ndim).astype(numpy.float64)
    B = numpy.random.randn(size*ndim).astype(numpy.float64)
    scale = numpy.float64(10.0)
    grad = numpy.zeros(size*ndim).astype(numpy.float64)
    cost = numpy.zeros(size).astype(numpy.float64)

    #time the reference function
    arguments = [cost, A, B, size, size, ndim, scale, grad]
    with open(get_kernel_path('gausstransform')+'gausstransform_c.cpp', 'r') as f:
        kernel_string = f.read()
    tune_params = {"block_size_x": [1]}
    #print("CPU timing")
    #tune_kernel("time_GaussTransform", kernel_string, size, arguments, tune_params,
    #            lang="C", compiler_options=['-I'+get_kernel_path('gausstransform'), '-O3'])

    #tune the GPU function
    print("GPU timing")

    with open(get_kernel_path('gausstransform')+'kernels.cu', 'r') as f:
        kernel_string = f.read()

    scale_sq = (scale*scale).astype(numpy.float64)
    arguments = [A, B, size, size, scale_sq, grad, cost]
    cp = ['-O3']

    print("GaussTransform kernel")
    tune_params = {"block_size_x": [32, 64, 128, 256, 512, 1024]}
    kernel1 = tune_kernel("GaussTransform", kernel_string, size, arguments, tune_params,
                grid_div_x=[], compiler_options=cp)

    print("Reduce cross term kernel")
    arguments = [numpy.zeros(1).astype(numpy.float64), cost, size, size, size]
    kernel2 = tune_kernel("reduce_cross_term", kernel_string, 1, arguments, tune_params,
                grid_div_x=[], compiler_options=cp)

    best_config1 = min(kernel1[0], key=lambda x:x['time'])
    best_config2 = min(kernel2[0], key=lambda x:x['time'])

    print("best GPU configuration, total time=", best_config1['time'] + best_config2['time'])
    print(best_config1)
    print(best_config2)



if __name__ == "__main__":
    tune_gausstransform()
