# This is a Makefile to build cln3del from Fortran sources using GNU Make.
# Ensure `make` is installed and in your shell's path by using `which`
# to find it's location. If the `make` command is not available,
# nothing will be output.
#
# > which make
# /usr/bin/make
#
# From the directory this script is in, run `make` and the first
# rule defined in the Makefile will be built.
# In this file, the `all` rule will build the `main` executable.

EXE = HILL_QNSTOP # Output executable.
SRC = model_hill.f95 main_HILL.f95 qnstops.f95 zigarray.f95 real_precision.f95
OBJ = $(SRC:.f95=.o) # Object files (source files with '.f95' replaced with '.o')

# Fortran 95 compiler and flags.
F95 = gfortran

# Add compiler specific flags.
ifeq ($(F95),g95)
	FFLAGS+=-O3 -Wall -ftrace=full -std=f95
else ifeq ($(F95),gfortran)
	FFLAGS+=-O3 -Wall
	FFLAGS+=-ffpe-trap=invalid,zero,overflow
	#FFLAGS+=-ffpe-summary=all
	FFLAGS+=-fcheck=all -fbacktrace
# Use the following line to enable OpenMP.
	FFLAGS+=-fopenmp
# For debugging, uncomment the following line.
	FFLAGS+=-g
else ifeq ($(F95),lf95)
	FFLAGS+=-Wa,--32 --f95 --in --trap
else ifeq ($(F95),f95)
	FFLAGS+=-en -O2 -s
	FFLAGS+=-TENV:simd_imask=off -TENV:simd_zmask=off -TENV:simd_omask=off
endif

# Libraries that need to be linked.
# LINK = -llapack -lblas 
# If using custom BLAS and LAPACK installations,
# uncomment the following line.
# LINK+=-L$(HOME)/BLAS -L$(HOME)/lapack-3.4.2

# case 1: works when LAPACK and BLAS are installed.
#OBJLIBS =
#LIBS    = -llapack -lblas

# case 2: works when LAPACK is not installed but BLAS is.
#LIBS    = -lblas
#OBJLIBS = lapack.o

# case 3: works when neither LAPACK nor BLAS are installed.
# LINK = blas.o lapack.o
# LIBS	=
SUP = blas_exp.o lapack_exp.o blas_qnstop.o lapack_qnstop.o expokit.o mataid.o

all: $(EXE)

$(EXE): $(OBJ) $(SUP)
	$(F95) $(FFLAGS) -o $@ $^
# real_precision.mod: real_precision.f95
# 	$(F95) $(FFLAGS) -c $< -o $(@:.mod=.o)
# # zigarray.mod: zigarray.f95 real_precision.mod
# # 	$(F95) $(FFLAGS) -c $< -o $(@:.mod=.o)
# qnstops.mod: qnstops.f95 zigarray.mod real_precision.mod
# 	$(F95) $(FFLAGS) -c $< -o $(@:.mod=.o)
# # model_hill.mod: model_hill.f95 real_precision.mod
# # 	$(F95) $(FFLAGS) -c $< -o $(@:.mod=.o)
# main_HILL.o: main_HILL.f95 model_hill.mod qnstops.mod 
# 	$(F95) $(FFLAGS) -c $< -o $@
# %.o %.mod: %.f95 real_precision.mod
# 	$(F95) $(FFLAGS) -c $< -o $(@:.mod=.o)



real_precision.o: real_precision.f95
	$(F95) $(FFLAGS) -c $< -o $@
zigarray.o: zigarray.f95 real_precision.o
	$(F95) $(FFLAGS) -c $< -o $@
qnstops.o: qnstops.f95 zigarray.o real_precision.o
	$(F95) $(FFLAGS) -c $< -o $@
model_hill.o: model_hill.f95 real_precision.o
	$(F95) $(FFLAGS) -c $< -o $@
main_HILL.o: main_HILL.f95 model_hill.o qnstops.o
	$(F95) $(FFLAGS) -c $< -o $@
%.o:%.f
	$(F95) -O3 -c $<


purge:	clean

clean:
		rm -f *~ *.o *.mod
	   	rm -f sOutput.txt

