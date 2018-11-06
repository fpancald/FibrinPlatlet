################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -static -std=c++11 -Wall -Wextra
NVCCFLAGS := 
IFLAGS = -isystem/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/
LIBS := -lgsl -lgslcblas -lpugixml -L/afs/crc.nd.edu/user/s/sbritto2/networkgen/FibrinCode/UnderConstructionThrustExtension/pugixml/lib64

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../PlateletForceDevice.cu \
../TorsionSolveOnDevice.cu \
../WLCSolveOnDevice.cu \
../AdvancePositionOnDevice.cu \
../LinkNodesOnDevice.cu \
../BucketSchemeOnDevice.cu \
../NodeSystemDevice.cu \
../NodeSystemBuilder.cpp \
../ForceDiagramStorage.cpp \
../main.cpp 


# this is a variable
OBJS += \
./PlateletForceDevice.o \
./TorsionSolveOnDevice.o \
./WLCSolveOnDevice.o \
./AdvancePositionOnDevice.o \
./LinkNodesOnDevice.o \
./BucketSchemeOnDevice.o \
./NodeSystemDevice.o \
./NodeSystemBuilder.o \
./ForceDiagramStorage.o \
./main.o 

 
CPP_DEPS += \
./PlateletForceDevice.d \
./TorsionSolveOnDevice.d \
./WLCSolveOnDevice.d \
./AdvancePositionOnDevice.d \
./LinkNodesOnDevice.d \
./BucketSchemeOnDevice.d \
./NodeSystemDevice.d \
./NodeSystemBuilder.d \
./ForceDiagramStorage.d \
./main.d 


#cpp files
%.o : ./%.cpp 
	$(CXX) $(CFLAGS) $(IFLAGS) -o  $@ -c $^

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $^ 



