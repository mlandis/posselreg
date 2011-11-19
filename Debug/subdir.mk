################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Alignment.cpp \
../ConTree.cpp \
../CondLikes.cpp \
../MbBitfield.cpp \
../MbEigensystem.cpp \
../MbMath.cpp \
../MbRandom.cpp \
../MbTransitionMatrix.cpp \
../Mcmc.cpp \
../Model.cpp \
../NodeKey.cpp \
../NodeSiteKey.cpp \
../NodeSiteVals.cpp \
../NodeVal.cpp \
../Parm.cpp \
../Parm_freqs.cpp \
../Parm_kappa.cpp \
../Parm_lambda.cpp \
../Parm_omega.cpp \
../Parm_omegawts.cpp \
../Parm_tree.cpp \
../Settings.cpp \
../TiProbs.cpp \
../TreeSum.cpp \
../Util.cpp \
../code.cpp \
../iomanager.cpp \
../main.cpp 

OBJS += \
./Alignment.o \
./ConTree.o \
./CondLikes.o \
./MbBitfield.o \
./MbEigensystem.o \
./MbMath.o \
./MbRandom.o \
./MbTransitionMatrix.o \
./Mcmc.o \
./Model.o \
./NodeKey.o \
./NodeSiteKey.o \
./NodeSiteVals.o \
./NodeVal.o \
./Parm.o \
./Parm_freqs.o \
./Parm_kappa.o \
./Parm_lambda.o \
./Parm_omega.o \
./Parm_omegawts.o \
./Parm_tree.o \
./Settings.o \
./TiProbs.o \
./TreeSum.o \
./Util.o \
./code.o \
./iomanager.o \
./main.o 

CPP_DEPS += \
./Alignment.d \
./ConTree.d \
./CondLikes.d \
./MbBitfield.d \
./MbEigensystem.d \
./MbMath.d \
./MbRandom.d \
./MbTransitionMatrix.d \
./Mcmc.d \
./Model.d \
./NodeKey.d \
./NodeSiteKey.d \
./NodeSiteVals.d \
./NodeVal.d \
./Parm.d \
./Parm_freqs.d \
./Parm_kappa.d \
./Parm_lambda.d \
./Parm_omega.d \
./Parm_omegawts.d \
./Parm_tree.d \
./Settings.d \
./TiProbs.d \
./TreeSum.d \
./Util.d \
./code.d \
./iomanager.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


