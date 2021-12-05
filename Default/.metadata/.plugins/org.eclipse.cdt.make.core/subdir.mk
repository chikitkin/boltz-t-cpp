################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../.metadata/.plugins/org.eclipse.cdt.make.core/specs.cpp 

C_SRCS += \
../.metadata/.plugins/org.eclipse.cdt.make.core/specs.c 

OBJS += \
./.metadata/.plugins/org.eclipse.cdt.make.core/specs.o 

CPP_DEPS += \
./.metadata/.plugins/org.eclipse.cdt.make.core/specs.d 

C_DEPS += \
./.metadata/.plugins/org.eclipse.cdt.make.core/specs.d 


# Each subdirectory must supply rules for building sources it contributes
.metadata/.plugins/org.eclipse.cdt.make.core/%.o: ../.metadata/.plugins/org.eclipse.cdt.make.core/%.c .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) oneAPI DPC++ Compiler'
	dpcpp -Wall -I/home/egor/git/boltz-t-cpp -I/opt/intel/oneapi/vpl/2021.4.0/include -I/opt/intel/oneapi/tbb/2021.3.0/env/../include -I/opt/intel/oneapi/mpi/2021.3.0//include -I/opt/intel/oneapi/mkl/2021.3.0/include -I/opt/intel/oneapi/ipp/2021.3.0/include -I/opt/intel/oneapi/ippcp/2021.3.0/include -I/opt/intel/oneapi/dpl/2021.4.0/linux/include -I/opt/intel/oneapi/dpcpp-ct/2021.3.0/include -I/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib -I/opt/intel/oneapi/dev-utilities/2021.3.0/include -I/opt/intel/oneapi/dal/2021.3.0/include -I/opt/intel/oneapi/compiler/2021.3.0/linux/include -I/opt/intel/oneapi/ccl/2021.3.0/include/cpu_gpu_dpcpp -I/opt/intel/oneapi/clck/2021.3.0/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/MEDCOUPLING/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/medfile/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/hdf5/include -I/home/egor/git/boltz-t-cpp -MMD -MP -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

.metadata/.plugins/org.eclipse.cdt.make.core/%.o: ../.metadata/.plugins/org.eclipse.cdt.make.core/%.cpp .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) oneAPI DPC++ Compiler'
	dpcpp -Wall -I/home/egor/git/boltz-t-cpp -I/opt/intel/oneapi/vpl/2021.4.0/include -I/opt/intel/oneapi/tbb/2021.3.0/env/../include -I/opt/intel/oneapi/mpi/2021.3.0//include -I/opt/intel/oneapi/mkl/2021.3.0/include -I/opt/intel/oneapi/ipp/2021.3.0/include -I/opt/intel/oneapi/ippcp/2021.3.0/include -I/opt/intel/oneapi/dpl/2021.4.0/linux/include -I/opt/intel/oneapi/dpcpp-ct/2021.3.0/include -I/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib -I/opt/intel/oneapi/dev-utilities/2021.3.0/include -I/opt/intel/oneapi/dal/2021.3.0/include -I/opt/intel/oneapi/compiler/2021.3.0/linux/include -I/opt/intel/oneapi/ccl/2021.3.0/include/cpu_gpu_dpcpp -I/opt/intel/oneapi/clck/2021.3.0/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/MEDCOUPLING/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/medfile/include -I/home/egor/MEDCOUPLING-9.7.0-MPI-UB20.04-SRC/BINARIES-UB20.04/hdf5/include -I/home/egor/git/boltz-t-cpp -MMD -MP -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


