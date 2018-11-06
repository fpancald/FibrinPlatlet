

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include subdir.mk
-include objects.mk


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs


# Add inputs and outputs from these tool invocations to the build variables

# All Target
all: bend-model

#Flags
CXXFLAGS=-O2 -std=c++0x -pg -g -c -Wall
NVCCFLAGS=-O2 -g -G
# Tool invocations
bend-model: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC C++ Linker'
	nvcc -o "bend-model" $(OBJS) $(USER_OBJS) $(ILIBS) $(LIBS) $(NVCCFLAGS)
	@echo 'Finished building target: $@'
	@echo ' '
#	$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) bend-model
	-@echo ' '

#post-build:
#	-mkdir --parents ../bend-0.0.6/Debug; cd ../../bend-0.0.6/Debug; cp ../../bend-model/Debug/bend-model bend-model
#	-@echo ' '

.PHONY: all clean dependents
.SECONDARY: post-build

-include ../makefile.targets
# DO NOT DELETE
