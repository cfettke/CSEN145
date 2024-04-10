#
# Instructions:
# You need to do three steps, see below
#

###################################################################
# Step 1: Configure your host compiler
###################################################################

# Serial/OpenMP
CC = g++
CFLAGS = -std=c++11 -O3 -march=native -mtune=native -funroll-loops -mfpmath=both

# OpenACC -- use -gpucc70 for V100 GPUs (gpu partition) and -gpu=cc61 for the GTX1080 GPUs (hub partition).
PGCC = pgc++
PGFLAGS = -O3 -acc -fast -gpu=cc61 -Minfo=accel

# CUDA
CUDA_CFLAGS = -std=c++11 -O3
CUDA_LDFLAGS =

# OpenMP
OMP_CFLAGS = -O3 -std=c++11 -foffload=-misa=sm_70 -fopenmp
OMP_LDFLAGS = 


RM = rm -f


###################################################################
# Step 2: Select which implementations to enable
################################################################### 
TARGETS = matmult_cpu \
matmult_omp_cpu \
matmult_omp_gpu \
matmult_acc \
matmult_hip \
matmult_cuda


all: $(TARGETS)

matmult_cpu: matmult_cpu.cpp matmult.h
	$(CC) $(CFLAGS) -o $(@) matmult_cpu.cpp

matmult_cuda: matmult_cuda.cu matmult.h
	nvc++ $(CUDA_CFLAGS) $(CUDA_LDFLAGS) -o $(@) matmult_cuda.cu
                
matmult_hip: matmult_hip.cpp matmult.h
	hipcc $(HIP_CFLAGS) $(HIP_LDFLAGS) -o $(@) matmult_hip.cpp
                
matmult_omp_cpu: matmult_omp_cpu.cpp matmult.h
	$(CC) $(OMP_CFLAGS) $(OMP_LDFLAGS) matmult_omp_cpu.cpp -o $(@) 
                
matmult_omp_gpu: matmult_omp_gpu.cpp matmult.h
	$(CC) $(OMP_CFLAGS) $(OMP_LDFLAGS) matmult_omp_gpu.cpp -o $(@) 

matmult_acc: matmult_acc.cpp matmult.h
	$(PGCC) $(PGFLAGS) -o $(@) matmult_acc.cpp
        
clean:
	$(RM) $(TARGETS) *.lib *.a *.exe *.obj *.o *.exp *.pyc

