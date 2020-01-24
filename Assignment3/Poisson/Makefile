TARGET_gpu1	= jacobi_gpu1
TARGET_gpu2	= jacobi_gpu2
TARGET_gpu3	= jacobi_gpu3
TARGET_gpu4	= jacobi_gpu4

LIBSRCS	= main.cu alloc3d_gpu.cu transfer3d_gpu.cu alloc3d.cu print.cu
LIBOBJS	= alloc3d_gpu.o transfer3d_gpu.o alloc3d.o print.o

MAIN_gpu1 = main_gpu1.o
MAIN_gpu2 = main_gpu2.o
MAIN_gpu3 = main_gpu3.o
MAIN_gpu4 = main_gpu4.o

OBSJ_gpu1 = $(MAIN_gpu1) jacobi_gpu1.o
OBSJ_gpu2 = $(MAIN_gpu2) jacobi_gpu2.o
OBSJ_gpu3 = $(MAIN_gpu3) jacobi_gpu3.o 
OBSJ_gpu4 = $(MAIN_gpu4) jacobi_gpu4.o

OPT	= -g -O3
PIC =
OMP   = -fopenmp
XPIC  = -Xcompiler
XOPT  = -Xptxas=-v -G # -lineinfo for profiler, use -G for debugging
XARCH = -arch=sm_70

CXX	= nvcc
CXXFLAGS = --compiler-options "$(OPT) $(PIC) $(OMP)" $(XARCH) $(XOPT) $(XPIC)

CUDA_PATH ?= /appl/cuda/10.2
INCLUDES = -I$(CUDA_PATH)/include \
	       -I$(CUDA_PATH)/samples/NVIDIA_CUDA-10.2_Samples/common/inc

SOFLAGS =
XLIBS	= -lcublas

all: $(TARGET_gpu1) $(TARGET_gpu2) $(TARGET_gpu3) $(TARGET_gpu4)

$(TARGET_gpu1): export JACOBIFLAG=-D_JACOBIGPU1
$(TARGET_gpu1): $(LIBOBJS) $(OBSJ_gpu1)
	$(CXX) -o $@  $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $^ $(XLIBS)

$(TARGET_gpu2): export JACOBIFLAG=-D_JACOBIGPU2
$(TARGET_gpu2): $(LIBOBJS) $(OBSJ_gpu2)
	$(CXX) -o $@ $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $^ $(XLIBS)

$(TARGET_gpu3): export JACOBIFLAG=-D_JACOBIGPU3
$(TARGET_gpu3): $(LIBOBJS) $(OBSJ_gpu3)
	$(CXX) -o $@ $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $^ $(XLIBS)

$(TARGET_gpu4): export JACOBIFLAG=-D_JACOBIGPU4
$(TARGET_gpu4): $(LIBOBJS) $(OBSJ_gpu4)
	$(CXX) -o $@ $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $^ $(XLIBS)



#all: $(MAIN_gpu1) $(MAIN_gpu2)

$(MAIN_gpu1):
	$(CXX) -o $@ -D_JACOBIGPU1 $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $(JACOBIFLAG) -c main.cu

$(MAIN_gpu2):
	$(CXX) -o $@ -D_JACOBIGPU2 $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $(JACOBIFLAG) -c main.cu

$(MAIN_gpu3):
	$(CXX) -o $@ -D_JACOBIGPU3 $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $(JACOBIFLAG) -c main.cu

$(MAIN_gpu4):
	$(CXX) -o $@ -D_JACOBIGPU4 $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $(JACOBIFLAG) -c main.cu

.SUFFIXES: .cu
.cu.o:
	$(CXX) -o $*.o -c $*.cu $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $(JACOBIFLAG) 

clean:
	/bin/rm -f $(TARGET_gpu1) $(TARGET_gpu2) $(LIBOBJS)
 
main_gpu1.o: main.cu alloc3d_gpu.h jacobi_gpu1.h transfer3d_gpu.h print.h
main_gpu2.o: main.cu alloc3d_gpu.h jacobi_gpu2.h transfer3d_gpu.h print.h
main_gpu3.o: main.cu alloc3d_gpu.h jacobi_gpu3.h transfer3d_gpu.h print.h
main_gpu4.o: main.cu alloc3d_gpu.h jacobi_gpu4.h transfer3d_gpu.h print.h
