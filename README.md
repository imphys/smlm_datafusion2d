# Software for Template-Free 2D Particle Fusion in Localization Microscopy.

This software implements a template-free particle fusion algorithm based on 
an all-to-all registration, which provides robustness against individual 
mis-registrations and underlabeling. The method does not assume any prior
knowledge about the structure to be reconstructed (template-free) and directly
works on localization data not pixelated images.

The main code is written in MATLAB and some of the compute-intensive kernels 
have been written in CUDA and C++.

## requirements
This code is built for a Linux enviroment. It might or might not work on a mac. 
For Windows a Linux shell enviroment could work.
For the CPU only code, no special libaries are needed for the compilation of 
the C code other than a C compiler.
For the GPU code, a CUDA compiler and libraries must be present and the CUB library
on the CPATH enviroment variable (see below).


## Installation

A Makefile is provided that can be used to compile the code and produce the 
required mex files.

Simply type ``make`` in the top-level directory and mex files will be 
produced in the `MATLAB/all2all/` directory.

If you do not have a CUDA-capable GPU use ``make cpu`` instead, to only
compile the CPU code.


## Example Usage
The DIPImage toolbox for MATLAB is required, please see http://www.diplib.org 
for installation instructions.

An example of how to use the code on experimental and simulated data is shown
in the MATLAB script `demo_all2all.m`. 


## Installation instructions for GPU Version

The mex files that call GPU functions will only be compiled by the Makefile 
if you have *nvcc* (the Nvidia CUDA compiler) installed. Just type `make` in 
the top-level directory after you've made sure that you've installed CUDA 
and the CUB library, see the instructions below.

### CUDA

The GPU code requires a CUDA-capable GPU as well as the CUDA toolkit to be 
installed. Please see Nvidia's website for installation fo the CUDA toolkit 
(https://developer.nvidia.com/cuda-downloads).

### CUB Library

The GPU code currently has one dependency, which is the CUB library. You can 
download it from: https://nvlabs.github.io/cub/index.html The easiest way to 
install CUB is to add the directory where you unpack CUB to your ``$CPATH`` 
environment variable.


## Troubleshooting

### Matlab mex headers not found

The Makefile tries to find the directories in which MATLAB was installed on 
your system. If this fails, you can manually insert the path to your MATLAB 
installation (ending with `/extern/include`) inside the Makefile. 

### CUDA headers not found

The Makefile also tries to automatically find the directories with headers 
and libraries needed to compile the CUDA codes. If this fails, these can as well be 
inserted at the top of the Makefile.

### <cub/cub.cuh> not found

The GPU code has only one external dependency, which is the CUB library. You 
can download it from: https://nvlabs.github.io/cub/index.html. The easiest 
way to install CUB is to add the top-level directory of where you've 
unpacked the CUB source codes to your ``$CPATH`` environment variable. For 
example, if you've unzipped the CUB sources into a directory called 
``/home/username/cub-version.number``, you can use 
``export CPATH=$CPATH:/home/username/cub-version.number/:`` to install CUB. In this way the 
nvcc compiler is able to find the CUB headers.

### Program tries to run GPU code when no GPU is present

Note that the mex files for the GPU code will be produced by `make` if your 
machine has `nvcc`. Once the mex files for the GPU code have been produced, 
the MATLAB code will prefer to use the GPU functions instead of the CPU 
functions. If you have no GPU available but did compile the mex files for 
the GPU code, you will get errors and MATLAB will exit. To disable the use 
of the GPU code type `make clean` and use `make cpu` instead of `make` or 
`make all`.

### Further questions

For further questions you can contact the authors:  
[Hamidreza Heydarian](https://hrheydarian.github.io/) (<H.Heydarian@tudelft.nl>) and  
[Ben van Werkhoven](https://github.com/benvanwerkhoven) (<b.vanwerkhoven@esciencecenter.nl>)  
[Bernd Rieger](http://homepage.tudelft.nl/z63s8/) (<b.rieger@tudelft.nl>)

Note that some files have been reused and adapted from the following sources:
GMM registration:  
    https://github.com/bing-jian/gmmreg  
    	[1] Jian, B. & Vemuri, B. C. Robust point set registration using Gaussian mixture models. IEEE PAMI 33, 16331645 (2011).  

Lie-algebraic averaging:  
    http://www.ee.iisc.ac.in/labs/cvl/research/rotaveraging/  
    [2] Govindu, V. Lie-algebraic averaging for globally consistent motion estimation. In Proc. IEEE Conf. on Computer Vision and Pattern Recognition (2004).  
    [3] Chatterjee, A. Geometric calibration and shape refinement for 3D reconstruction PhD thesis. Indian Institute of Science (2015).
l1-magic optimization toolbox:  
    https://statweb.stanford.edu/~candes/l1magic/
Natural-Order Filename Sort  
    https://nl.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

## Developer instructions

The testing and tuning scripts for the GPU code have been written in Python, 
using [Kernel Tuner](https://github.com/benvanwerkhoven/kernel_tuner). This 
section provides information on how to setup a development environment. Note 
that these steps are only needed if you are interested in modifying the CUDA 
and C++ codes.

### Python 3

The tests for the GPU code and several of the C functions are written in 
Python, to run these a Python 3 installation is required. The easiest way to 
get this is using [Miniconda](https://conda.io/miniconda.html).

On Linux systems one could type the following commands to download and 
install Python 3 using Miniconda:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

All the required Python packages can be installed using the following command,
before you run it make sure that you have CUDA installed:
```
pip install -r requirements.txt
```

The tests can be run using ``nose``, for example by typing the following in 
the top-level or test directory:
```
nosetests -v
```
