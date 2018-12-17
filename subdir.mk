################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -static -std=c++11 -Wall -Wextra
NVCCFLAGS :=
IFLAGS = -isystem/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/

CURRENT_DIR := $(shell pwd)
LIBS:=  -lgsl -lgslcblas -lpugixml -L/$(CURRENT_DIR)/pugixml/lib64

# Add inputs and outputs from these tool invocations to the build variables
CPP_SRCS += \
../Plt_Field_Node_Force.cu \
../Plt_Field_Plt_Force.cu \
../Plt_Arm_Node_Force.cu \
../Plt_Arm_Plt_Force.cu \
../Plt_Vol_Exc_Force.cu \
../Torsion_Force.cu \
../WLC_Force.cu \
../Advance_Positions.cu \
../Link_Nodes.cu \
../Bucket_Net.cu \
../Bucket_Plt.cu \
../System.cu \
../System_Builder.cpp \
../Storage.cpp \
../main.cpp


# this is a variable
OBJS += \
./Plt_Field_Node_Force.o \
./Plt_Field_Plt_Force.o \
./Plt_Arm_Node_Force.o \
./Plt_Arm_Plt_Force.o \
./Plt_Vol_Exc_Force.o \
./Torsion_Force.o \
./WLC_Force.o \
./Advance_Positions.o \
./Link_Nodes.o \
./Bucket_Net.o \
./Bucket_Plt.o \
./System.o \
./System_Builder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./Plt_Field_Node_Force.d \
./Plt_Field_Plt_Force.d \
./Plt_Arm_Node_Force.d \
./Plt_Arm_Plt_Force.d \
./Plt_Vol_Exc_Force.d \
./Torsion_Force.d \
./WLC_Force.d \
./Advance_Positions.d \
./Link_Nodes.d \
./Bucket_Net.d \
./Bucket_Plt.d \
./System.d \
./System_Builder.d \
./Storage.d \
./main.d


#cpp files
%.o : ./%.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) -o  $@ -c $^


#cuda files
%.o : ./%.cu
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $^




