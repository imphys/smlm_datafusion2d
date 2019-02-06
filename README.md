# Software for Template-Free 2D Particle Fusion in Localization Microscopy.

This software implements a template-free particle fusion algorithm based on 
an all-to-all registration, which provides robustness against individual 
mis-registrations and underlabeling. The method does not assume any prior
knowledge about the structure to be reconstructed (template-free) and directly
works on localization data not pixelated images.

The main code is written in MATLAB and some of the compute-intensive kernels 
have been written in CUDA and C++.

## Installation on Linux

Use the following commands to build the necessary libraries for this software:

```bash

cmake .
make
make install
````

CMake locates the code's dependencies and generates a Makefile. Make compiles the mex files and necessary shared libraries. Make install copies the mex files to the right directory for use in Matlab, it does not require priviledges to run.

## Installation on Windows

Make sure you have all dependencies downloaded and installed. In particular make sure that your CUDA version and MS Visual Studio version are 
compatible. Use the following commands to build the necessary libraries for this software:

```bash

mkdir build
cd build
cmake -G "Visual Studio 14 2015 Win64" -DCUB_INCLUDE_DIR="C:\path_to_your_cub_dir\cub" ..
cd ..
python fix_msbuild_files.py

````

Then open MS Visual Studio, for example by double clicking on ALL_BUILD.vcxproj and build all the targets.
It could be that you manually need to build the targets expdist and gausstransform before building the 
other targets.

Then, copy all the .mexw64 files from build/Debug/ to the MATLAB/all2all directory.



### Installation instructions for CPU-only Version

If you do not have a CUDA capable GPU you could install the software without GPU acceleration. Do note that the code will be orders of magnitude slower. Use the following commands to install the CPU-only version:
```bash

cmake -DUSE_GPU=OFF .
make
make install
```

## Example Usage

An example of how to use the code on experimental and simulated data is shown
in the MATLAB script `demo_all2all.m`. 

## Dependencies

### CMake

CMake is used as the build system for compiling the code. Here is a quick way to install the latest version of CMake locally without the need for root access:  
1. Run the following from your terminal and also add it to your ~/.bashrc file for future usage. 
```bash
export PATH=$HOME/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lib/:$LD_LIBRARY_PATH
```
2. Head over to the [CMake downloads page](https://cmake.org/download/) and get the latest “Unix/Linux Source” \*.tar.gz file.  
3. Run the following command sequentially (this can take a few minutes !):  
```bash

tar -xf cmake*.tar.gz
cd cmake*
./configure --prefix=$HOME
make
make install
```  
4. You should now have the new installation of cmake ready. Check the version by:  
```bash
cmake --version
```

### MATLAB

The main code is written in Matlab, it's not possible to use the software without it.
The DIPImage toolbox for MATLAB is required, please see http://www.diplib.org for installation instructions.

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

### Operation System

This code has been developed for a Linux enviroment.
We have yet to test it on OSX and Windows, but the CMake build system should support compilation on these platforms.
If you run into issues using the software on Windows or Mac please create an issue on GitHub or contact the authors.

### CUB library or <cub/cub.cuh> not found

The GPU code has only one external dependency, which is the CUB library. 
cmake will try to find and locate the CUB library on your system. If that fails it will
attempt to download the library directly. If that fails, you could install the library manually.

You can download CUB from: https://nvlabs.github.io/cub/index.html. The easiest 
way to install CUB is to add the top-level directory of where you've 
unpacked the CUB source codes to your ``$CPATH`` environment variable. For 
example, if you've unzipped the CUB sources into a directory called 
``/home/username/cub-version.number``, you can use 
``export CPATH=$CPATH:/home/username/cub-version.number/:`` to install CUB. In this way the 
nvcc compiler is able to find the CUB headers.

### Program tries to run GPU code when no GPU is present

Note that the mex files for the GPU code will be installed into the MATLAB/all2all directory.  
Once the mex files for the GPU code have been produced, 
the MATLAB code will prefer to use the GPU functions instead of the CPU 
functions. If you have no GPU available but did compile the mex files for 
the GPU code, you will get errors and MATLAB will exit. To disable the use 
of the GPU code remove the mex files from the MATLAB/all2all directory and reinstall using ``cmake -DUSE_GPU=OFF .``.

### Further questions

For further questions feel free to create an issue on GitHub. You can also contact the authors:  
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
## Reference

If you find this code useful for your research, please cite

```
@article {PMID:30224671,
	Title = {Template-free 2D particle fusion in localization microscopy},
	Author = {Heydarian, Hamidreza and Schueder, Florian and Strauss, Maximilian T and van Werkhoven, Ben and Fazel, Mohamadreza and Lidke, Keith A and Jungmann, Ralf and Stallinga, Sjoerd and Rieger, Bernd},
	DOI = {10.1038/s41592-018-0136-6},
	Number = {10},
	Volume = {15},
	Month = {October},
	Year = {2018},
	Journal = {Nature methods},
	ISSN = {1548-7091},
	Pages = {781—784},
	URL = {https://doi.org/10.1038/s41592-018-0136-6},
}
```
