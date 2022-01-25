################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cli.cpp \
../src/full.cpp \
../src/mesh-test.cpp \
../src/mesh.cpp \
../src/mesh2.cpp \
../src/run.cpp \
../src/solver.cpp \
../src/tensor-test.cpp \
../src/tucker.cpp 

OBJS += \
./src/cli.o \
./src/full.o \
./src/mesh-test.o \
./src/mesh.o \
./src/mesh2.o \
./src/run.o \
./src/solver.o \
./src/tensor-test.o \
./src/tucker.o 

CPP_DEPS += \
./src/cli.d \
./src/full.d \
./src/mesh-test.d \
./src/mesh.d \
./src/mesh2.d \
./src/run.d \
./src/solver.d \
./src/tensor-test.d \
./src/tucker.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) oneAPI DPC++ Compiler'
	dpcpp -Wall -I/home/egor/git/boltz-t-cpp -I/opt/intel/oneapi/vpl/2021.4.0/include -I/opt/intel/oneapi/tbb/2021.3.0/env/../include -I/opt/intel/oneapi/mpi/2021.3.0//include -I/opt/intel/oneapi/mkl/2021.3.0/include -I/opt/intel/oneapi/ipp/2021.3.0/include -I/opt/intel/oneapi/ippcp/2021.3.0/include -I/opt/intel/oneapi/dpl/2021.4.0/linux/include -I/opt/intel/oneapi/dpcpp-ct/2021.3.0/include -I/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib -I/opt/intel/oneapi/dev-utilities/2021.3.0/include -I/opt/intel/oneapi/dal/2021.3.0/include -I/opt/intel/oneapi/compiler/2021.3.0/linux/include -I/opt/intel/oneapi/ccl/2021.3.0/include/cpu_gpu_dpcpp -I/opt/intel/oneapi/clck/2021.3.0/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/MEDCOUPLING/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/medfile/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/hdf5/include -I/home/egor/git/boltz-t-cpp -MMD -MP -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


