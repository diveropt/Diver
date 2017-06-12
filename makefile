# ===========================================
# Makefile for Diver
# Author: Pat Scott, p.scott@imperial.ac.uk
#
# The Diver make system is not really complex enough to require autotools or
# cmake.  Just change this makefile by hand to suit your system, or call it
# from another makefile or the command line with overrides.
#
# The variables to override if you are just making
# libdiver are FF, FOPT and MODULE.
# If you also want to make the C/C++ examples, you may
# also need to override CC, COPT, MIXOPT_C and MIXOPT_CPP.
#
# Some typical examples are
#  Without MPI:
#    GNU:   make FF=gfortran CC=gcc MODULE=J FOPT="-ffixed-line-length-none -cpp"
#    intel: make FF=ifort CC=icc MODULE=module FOPT=-extend_source MIXOPT_C=-nofor-main MIXOPT_C=-nofor-main
#  With MPI:
#    GNU:   make FF=mpif90 CC=mpicc MODULE=J FOPT="-ffixed-line-length-none -DMPI -cpp" COPT=-DMPI
#    intel: make FF=mpif90 CC=mpicc MODULE=module FOPT="-extend_source -DMPI -fpp" COPT=-DMPI MIXOPT_C=-nofor-main MIXOPT_C=-nofor-main
#
#============================================

# Define fortran compiler and options
FF=mpif90
# Intel
#FOPT=-extend_source -DMPI -fpp #-warn all -check all -check noarg-temp_created
#MODULE=module
# GNU
FOPT=-ffixed-line-length-none -DMPI -cpp #-Wall -fcheck=all
MODULE=J

# Define C/C++ compilers and options
CC=mpicc
COPT=-DMPI
CPPOPT=
SO_LINK_FLAGS=
# Intel
#MIXOPT_C=-nofor-main
#MIXOPT_CPP=-nofor-main
# GNU
MIXOPT_C=
MIXOPT_CPP=

# Define some more internal variables
CLEAN = _clean
AR = ar r
SHAR = ld -shared
LIBNAME = diver
EXAMPLENAMES = example_f example_c example_cpp
EXAMPLECLEAN = $(EXAMPLENAMES:%=%$(CLEAN))
PREFIX = diver
DIVER_ROOT = $(PWD)
SOURCE = ${DIVER_ROOT}/src
INC = ${DIVER_ROOT}/include
BUILD = ${DIVER_ROOT}/build
LIB = ${DIVER_ROOT}/lib

# Set internal compile commands
DIVER_FF=$(FF)
DIVER_FOPT=$(FOPT) -O3 -fPIC -I$(INC) -$(MODULE) $(BUILD)
DIVER_CC=$(CC)
DIVER_COPT=$(COPT) -O3 -fPIC -I$(INC)
DIVER_CPPOPT=$(CPPOPT) -lstdc++
DIVER_MIXOPT_C=$(MIXOPT_C)
DIVER_MIXOPT_CPP=$(MIXOPT_CPP) -lstdc++
DIVER_SO_LINK_FLAGS=$(SO_LINK_FLAGS) -shared

# Send compile commands over to example makefiles
export DIVER_FF DIVER_FOPT DIVER_CC DIVER_COPT DIVER_CPPOPT DIVER_MIXOPT_C DIVER_MIXOPT_CPP LIBNAME PREFIX LIB INC

# Identify source files
SOURCEFILES = detypes deutils mutation crossover selection init converge posterior evidence io de cwrapper
OBJ = $(SOURCEFILES:%=$(BUILD)/%.o)

all: libdiver.a $(EXAMPLENAMES)

libdiver.a: makefile $(OBJ)
	$(AR) $(LIB)/$@ $(OBJ)

libdiver.so: makefile $(OBJ)
	$(DIVER_FF) $(DIVER_SO_LINK_FLAGS) -o $(LIB)/$@ $(TYPEOBJ) $(OBJ)

$(BUILD)/converge.o: $(SOURCE)/converge.f90 $(BUILD)/detypes.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/crossover.o: $(SOURCE)/crossover.f90 $(BUILD)/detypes.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/cwrapper.o: $(SOURCE)/cwrapper.f90 $(BUILD)/detypes.o $(BUILD)/de.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/de.o: $(SOURCE)/de.f90 $(BUILD)/detypes.o $(BUILD)/deutils.o $(BUILD)/init.o $(BUILD)/io.o $(BUILD)/converge.o $(BUILD)/selection.o $(BUILD)/mutation.o $(BUILD)/crossover.o $(BUILD)/posterior.o $(BUILD)/evidence.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/detypes.o: $(SOURCE)/detypes.f90
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/deutils.o: $(SOURCE)/deutils.f90 $(BUILD)/detypes.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/evidence.o: $(SOURCE)/evidence.f90 $(BUILD)/detypes.o $(BUILD)/deutils.o $(BUILD)/posterior.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/init.o: $(SOURCE)/init.f90 $(BUILD)/detypes.o $(BUILD)/deutils.o $(BUILD)/mutation.o $(BUILD)/crossover.o $(BUILD)/selection.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/io.o: $(SOURCE)/io.f90 $(BUILD)/detypes.o $(BUILD)/deutils.o $(BUILD)/evidence.o $(BUILD)/converge.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/mutation.o: $(SOURCE)/mutation.f90 $(BUILD)/detypes.o $(BUILD)/deutils.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/posterior.o: $(SOURCE)/posterior.f90 $(BUILD)/detypes.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(BUILD)/selection.o: $(SOURCE)/selection.f90 $(BUILD)/detypes.o $(BUILD)/deutils.o $(BUILD)/mutation.o $(BUILD)/crossover.o
	$(DIVER_FF) $(DIVER_FOPT) -c $< -o $@

$(EXAMPLENAMES): libdiver.a
	$(MAKE) EXAMPLENAME=$@ -C $@

clean:
	rm -f $(LIB)/*.a $(LIB)/*.so; \
	cd $(BUILD); rm -f *.o *.mod; \
	cd $(INC); rm -f *.mod

$(EXAMPLECLEAN):
	make EXAMPLENAME=$(subst $(CLEAN),,$@) \
         -C $(subst $(CLEAN),,$@) clean

cleanall: clean $(EXAMPLECLEAN)
