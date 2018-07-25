
# $< first dependency in the list
# $@ name of the target
# $^ all dependencies of this target
.PHONY: clean all cpu

####################
#
# You can manually insert the required paths for your setup:

# The MATLAB extern/include path (directory contains mex.h)
MATLAB_INC_PATH = /opt/ud/Matlab-R2016b/extern/include/

# The CUDA headers include path (directory contains cuda.h)
CUDA_INC_PATH = /opt/ud/cuda-8.0/include/

# The CUDA library include path (directory contains libcudart.so)
CUDA_LIB_PATH = /opt/ud/cuda-8.0/lib64/

# cub headers
# CUB_INC = /home/hheydarian/cub/
####################
#
# Do not change anything below this line

# Attempt to auto-detect the Matlab include path if none is specified above
ifeq ($(MATLAB_INC_PATH),)
MATLAB_RESULT := $(shell which matlab 2> NULL)
MATLAB_INC_PATH = $(dir $(MATLAB_RESULT))../extern/include
endif

# Attempt to auto-detect the CUDA include paths if none are specified above
NVCC_RESULT := $(shell which nvcc 2> NULL)
ifeq ($(CUDA_INC_PATH),)
CUDA_INC_PATH = $(dir $(NVCC_RESULT))../include/
endif
ifeq ($(CUDA_LIB_PATH),)
CUDA_LIB_PATH = $(dir $(NVCC_RESULT))../lib64/
endif

CPATH+=:$(CUDA_INC_PATH)

NVCC_FLAGS = -O3 -m64 -Xcompiler=-fPIC -Xptxas=-v
CC_FLAGS = -O3 -m64 -fPIC

SMS ?= 30 35 37 50 52 60
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))

LIBDIR = lib/

TARGET_DIR = MATLAB/all2all/

TARGETS_GPU = 
TARGETS_CPU = $(TARGET_DIR)/mex_gausstransform_cpu.mexa64 $(TARGET_DIR)/mex_expdist_cpu.mexa64

# Check if nvcc is installed and if so add the GPU targets
NVCC_TEST := $(notdir $(NVCC_RESULT))
ifeq ($(NVCC_TEST),nvcc)
TARGETS_GPU += $(TARGET_DIR)/mex_gausstransform.mexa64 $(TARGET_DIR)/mex_expdist.mexa64
endif

TARGETS = $(TARGETS_CPU) $(TARGETS_GPU)

export


all: $(TARGETS)

cpu: $(TARGETS_CPU)


$(TARGET_DIR)/mex_gausstransform.mexa64: gausstransform/mex_gausstransform.cpp $(LIBDIR)/gausstransform.so
	mex $^ -Igausstransform/src/ -I$(CUDA_INC_PATH) -L$(CUDA_LIB_PATH) -lcudart -output $@

$(TARGET_DIR)/mex_gausstransform_cpu.mexa64: gausstransform/mex_gausstransform_cpu.cpp
	mex $^ -output $@

$(TARGET_DIR)/mex_expdist_cpu.mexa64: expdist/mex_expdist_cpu.cpp
	mex $^ -output $@

$(LIBDIR)/gausstransform.so: gausstransform/gausstransform.cu gausstransform/kernels.cu gausstransform/gausstransform.h
	mkdir -p $(LIBDIR)
	nvcc $(NVCC_FLAGS) $(GENCODE_FLAGS) -c $< -I$(MATLAB_INC_PATH) -shared -o $@

$(TARGET_DIR)/mex_expdist.mexa64: expdist/mex_expdist.cpp $(LIBDIR)/expdist.so
	mex $^ -Iexpdist/src/ -I$(CUDA_INC_PATH) -L$(CUDA_LIB_PATH) -lcudart -output $@

$(LIBDIR)/expdist.so: expdist/expdist.cu expdist/kernels.cu expdist/expdist.h
	mkdir -p $(LIBDIR)
	nvcc $(NVCC_FLAGS) $(GENCODE_FLAGS) -c $< -I$(MATLAB_INC_PATH) -shared -o $@



clean:
	rm -rf $(TARGET_DIR)/*.mexa64 $(LIBDIR)
