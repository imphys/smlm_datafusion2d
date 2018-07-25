#!/usr/bin/env python

from collections import OrderedDict
import json
import numpy

from kernel_tuner import tune_kernel

from tune_utils import get_kernel_path


def tune_expdist():

    device = 2

    tune_params = OrderedDict()
    tune_params["block_size_x"] = [32] #[2**i for i in range(5,10)]
    tune_params["block_size_y"] = [2**i for i in range(6)]
    tune_params["tile_size_x"] = [2**i for i in range(4)]
    tune_params["tile_size_y"] = [2**i for i in range(4)]
    tune_params["use_shared_mem"] = [1] #[0, 1]

    #setup test input
    alloc_size = 22000
    size = numpy.int32(20000)
    max_blocks = numpy.int32( numpy.ceil(size / float(numpy.amin(tune_params["block_size_x"]))) *
                              numpy.ceil(size / float(numpy.amin(tune_params["block_size_y"]))) )
    ndim = numpy.int32(2)

    A = numpy.random.randn(alloc_size*ndim).astype(numpy.float64)
    B = A+0.00001*numpy.random.randn(alloc_size*ndim).astype(numpy.float64)
    scale_A = numpy.absolute(0.01*numpy.random.randn(alloc_size).astype(numpy.float64))
    scale_B = numpy.absolute(0.01*numpy.random.randn(alloc_size).astype(numpy.float64))
    cost = numpy.zeros((max_blocks)).astype(numpy.float64)

    #time the reference function
    #arguments = [cost, A, B, size, size, ndim, scale_A, scale_B]
    #with open(get_kernel_path('expdist')+'expdist_c.cpp', 'r') as f:
    #    kernel_string = f.read()
    #print("CPU timing")
    #tune_kernel("time_expdist", kernel_string, size, arguments, {"block_size_x": [1]},
    #           lang="C", compiler_options=['-I'+get_kernel_path('expdist'), '-O3'], device=2)

    #tune the GPU function
    print("GPU timing")

    with open(get_kernel_path('expdist')+'kernels.cu', 'r') as f:
        kernel_string = f.read()

    arguments = [A, B, size, size, scale_A, scale_B, cost]
    cp = ['-O3']

    grid_div_x = ["block_size_x", "tile_size_x"]
    grid_div_y = ["block_size_y", "tile_size_y"]

    kernel1 = tune_kernel("ExpDist", kernel_string, (size, size), arguments, tune_params,
               compiler_options=cp, grid_div_x=grid_div_x, grid_div_y=grid_div_y, device=2)

    with open("expdist.json", 'w') as fp:
        json.dump(kernel1, fp)

    best_config1 = min(kernel1[0], key=lambda x:x['time'])

    nblocks = numpy.int32( numpy.ceil(size / float(best_config1["block_size_x"]*best_config1["tile_size_x"])) *
                           numpy.ceil(size / float(best_config1["block_size_y"]*best_config1["tile_size_y"])) )

    tune_params = OrderedDict()
    tune_params["block_size_x"] = [32*i for i in range(1,33)]

    arguments = [numpy.zeros(1).astype(numpy.float64), cost, size, size, nblocks]
    kernel2 = tune_kernel("reduce_cross_term", kernel_string, 1, arguments, tune_params,
                grid_div_x=[], compiler_options=cp, device=2)

    best_config2 = min(kernel2[0], key=lambda x:x['time'])

    print("best GPU configuration, total time=", best_config1['time'] + best_config2['time'])
    print(best_config1)
    print(best_config2)

if __name__ == "__main__":
    tune_expdist()
