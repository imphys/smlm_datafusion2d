#!/usr/bin/env python

from collections import OrderedDict
import json
import numpy

from kernel_tuner import tune_kernel

from tune_utils import get_kernel_path


def tune_expdist():

    tune_params = OrderedDict()
    tune_params["block_size_x"] = [2**i for i in range(5,10)]
    tune_params["block_size_y"] = [2**i for i in range(6)]
    tune_params["tile_size_x"] = [2**i for i in range(4)]
    tune_params["tile_size_y"] = [2**i for i in range(4)]
    tune_params["use_shared_mem"] = [1] #[0, 1]

    #setup test input
    alloc_size = 22000
    size = numpy.int32(20000)
    max_blocks = numpy.int32( numpy.ceil(size / float(numpy.amin(tune_params["block_size_x"]))) )

    ndim = numpy.int32(2)

    A = numpy.random.randn(alloc_size*ndim).astype(numpy.float64)
    B = A+0.00001*numpy.random.randn(alloc_size*ndim).astype(numpy.float64)
    scale_A = numpy.absolute(0.01*numpy.random.randn(alloc_size).astype(numpy.float64))
    scale_B = numpy.absolute(0.01*numpy.random.randn(alloc_size).astype(numpy.float64))
    cost = numpy.zeros((max_blocks)).astype(numpy.float64)

    #tune the GPU function
    with open(get_kernel_path('expdist')+'kernels.cu', 'r') as f:
        kernel_string = f.read()

    arguments = [A, B, size, size, scale_A, scale_B, cost]
    cp = ['-O3']

    grid_div_x = ["block_size_x", "tile_size_x"]

    kernel1, env = tune_kernel("ExpDist_column", kernel_string, size, arguments, tune_params,
               compiler_options=cp, grid_div_x=grid_div_x)

    devname = "".join(env["device_name"].split()) #device name without whitespace

    with open("expdist_column_" + devname + ".json", 'w') as fp:
        json.dump([kernel1, env], fp)

    best_config1 = min(kernel1, key=lambda x:x['time'])

    nblocks = numpy.int32( numpy.ceil(size / float(best_config1["block_size_x"]*best_config1["tile_size_x"])) )

    tune_params = OrderedDict()
    tune_params["block_size_x"] = [32*i for i in range(1,33)]

    arguments = [numpy.zeros(1).astype(numpy.float64), cost, size, size, nblocks]
    kernel2 = tune_kernel("reduce_cross_term", kernel_string, 1, arguments, tune_params,
                grid_div_x=[], compiler_options=cp)

    best_config2 = min(kernel2[0], key=lambda x:x['time'])

    print("best GPU configuration, total time=", best_config1['time'] + best_config2['time'])
    print(best_config1)
    print(best_config2)

if __name__ == "__main__":
    tune_expdist()
