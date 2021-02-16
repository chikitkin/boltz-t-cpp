################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/mesh.cpp \
../src/solver.cpp \
../src/tensor_class.cpp 

O_SRCS += \
../src/mesh.o \
../src/tensor_class.o 

OBJS += \
./src/mesh.o \
./src/solver.o \
./src/tensor_class.o 

CPP_DEPS += \
./src/mesh.d \
./src/solver.d \
./src/tensor_class.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) C++ Compiler Classic'
	icpc -g -O0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


