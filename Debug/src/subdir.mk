################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/solver.cpp \
../src/tensor_class.cpp 

OBJS += \
./src/solver.o \
./src/tensor_class.o 

CPP_DEPS += \
./src/solver.d \
./src/tensor_class.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -g -O3 -mkl=sequential -I/home/egor/intel/compilers_and_libraries_2019.0.117/linux/mkl/include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


