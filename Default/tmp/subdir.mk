################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tmp/driver.cpp \
../tmp/init.cpp \
../tmp/normalize.cpp \
../tmp/qr.cpp \
../tmp/sum.cpp \
../tmp/svd_trunc.cpp \
../tmp/test_mesh.cpp 

OBJS += \
./tmp/driver.o \
./tmp/init.o \
./tmp/normalize.o \
./tmp/qr.o \
./tmp/sum.o \
./tmp/svd_trunc.o \
./tmp/test_mesh.o 

CPP_DEPS += \
./tmp/driver.d \
./tmp/init.d \
./tmp/normalize.d \
./tmp/qr.d \
./tmp/sum.d \
./tmp/svd_trunc.d \
./tmp/test_mesh.d 


# Each subdirectory must supply rules for building sources it contributes
tmp/%.o: ../tmp/%.cpp tmp/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) oneAPI DPC++ Compiler'
	dpcpp -Wall -I/home/egor/git/boltz-t-cpp -I/opt/intel/oneapi/vpl/2021.4.0/include -I/opt/intel/oneapi/tbb/2021.3.0/env/../include -I/opt/intel/oneapi/mpi/2021.3.0//include -I/opt/intel/oneapi/mkl/2021.3.0/include -I/opt/intel/oneapi/ipp/2021.3.0/include -I/opt/intel/oneapi/ippcp/2021.3.0/include -I/opt/intel/oneapi/dpl/2021.4.0/linux/include -I/opt/intel/oneapi/dpcpp-ct/2021.3.0/include -I/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib -I/opt/intel/oneapi/dev-utilities/2021.3.0/include -I/opt/intel/oneapi/dal/2021.3.0/include -I/opt/intel/oneapi/compiler/2021.3.0/linux/include -I/opt/intel/oneapi/ccl/2021.3.0/include/cpu_gpu_dpcpp -I/opt/intel/oneapi/clck/2021.3.0/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/MEDCOUPLING/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/medfile/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/hdf5/include -I/home/egor/git/boltz-t-cpp -MMD -MP -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


