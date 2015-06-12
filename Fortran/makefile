# ===========================================
# Makefile for Diver 
# Author: Pat Scott, p.scott@imperial.ac.uk
#
# The Diver make system is not really complex enough to require autotools or 
# cmake.  Just change this makefile by hand to suit your system, or call it
# from another makefile or the command line with overrides.  
# 
# The variables to override if you are just making
# libdiver are DIVER_FF, DIVER_FOPT and DIVER_MODULE.
# If you also want to make the C/C++ examples, you may
# also need to override DIVER_CC, DIVER_COPT, DIVER_MIXOPT_C and DIVER_MIXOPT_CPP.
#
# Some typical examples are
#  Without MPI:
#    GNU:   make DIVER_FF=gfortran DIVER_CC=gcc DIVER_MODULE=J DIVER_FOPT="-ffree-line-length-none -ffixed-line-length-none -cpp"
#    intel: make DIVER_FF=ifort DIVER_CC=icc DIVER_MODULE=module DIVER_FOPT=-extend_source DIVER_MIXOPT_C=-nofor-main DIVER_MIXOPT_C=-nofor-main
#  With MPI:
#    GNU:   make DIVER_FF=mpif90 DIVER_CC=mpicc DIVER_MODULE=J DIVER_FOPT="-ffree-line-length-none -ffixed-line-length-none -cpp -DMPI" DIVER_COPT=-DMPI
#    intel: make DIVER_FF=mpif90 DIVER_CC=mpicc DIVER_MODULE=module DIVER_FOPT="-extend_source -DMPI -fpp" DIVER_COPT=-DMPI DIVER_MIXOPT_C=-nofor-main DIVER_MIXOPT_C=-nofor-main
#
#============================================

# Define fortran compiler and options
FF=mpif90
FOPT=-extend_source -DMPI -warn all -check all -check noarg-temp_created -fpp
MODULE=module

# Define C/C++ compilers and options
CC=mpicc
COPT=-DMPI
CPPOPT=
MIXOPT_C=-nofor-main
MIXOPT_CPP=-nofor-main

# Define some more internal variables
CLEAN = _clean
AR = ar r  
SHAR = ld -shared
LIBNAME = diver
EXAMPLENAMES = example_f example_c example_cpp
EXAMPLECLEAN = $(EXAMPLENAMES:%=%$(CLEAN))
PREFIX = diver
DIVER_ROOT = $(PWD)
SOURCE = ${DIVER_ROOT}/source
INC = ${DIVER_ROOT}/include
BUILD = ${DIVER_ROOT}/build
LIB = ${DIVER_ROOT}/lib

# Set internal compile commands
FF=$(DIVER_FF)
FOPT=$(DIVER_FOPT) -O3 -fPIC -I$(INC) -$(DIVER_MODULE) $(BUILD)
CC=$(DIVER_CC)
COPT=$(DIVER_COPT) -O3 -fPIC -I$(INC)
CPPOPT=$(DIVER_CPPOPT) -lstdc++
MIXOPT_C=$(DIVER_MIXOPT_C)
MIXOPT_CPP=$(DIVER_MIXOPT_CPP) -lstdc++

# Send compile commands over to example makefiles
export FF FOPT CC COPT CPPOPT MIXOPT_C MIXOPT_CPP LIBNAME PREFIX LIB INC

# Identify source files
TYPESOURCE = detypes
TYPEOBJ = $(TYPESOURCE:%=$(BUILD)/%.o)
TYPEF90 = $(TYPESOURCE:%=$(SOURCE)/%.f90)
SOURCEFILES = deutils mutation crossover selection init converge posterior evidence io de cwrapper
OBJ = $(SOURCEFILES:%=$(BUILD)/%.o)

all: libdiver.a $(EXAMPLENAMES)

libdiver.a: makefile $(TYPEOBJ) $(OBJ)
	$(AR) $(LIB)/$@ $(TYPEOBJ) $(OBJ)

libdiver.so: makefile $(TYPEOBJ) $(OBJ) 
	$(SHAR) -o $(LIB)/$@ $(TYPEOBJ) $(OBJ)  
 
$(BUILD)/%.o: $(SOURCE)/%.f90 $(TYPEF90)
	$(FF) $(FOPT) -c $< -o $@

$(BUILD)/%.o: $(SOURCE)/%.F90 $(TYPEF90)
	$(FF) $(FOPT) -c $< -o $@

$(EXAMPLENAMES): libdiver.a
	make EXAMPLENAME=$@ -C $@

clean:
	rm -f $(LIB)/*.a $(LIB)/*.so; \
	cd $(BUILD); rm -f *.o *.mod; \
	cd $(INC); rm -f *.mod

$(EXAMPLECLEAN):
	make EXAMPLENAME=$(subst $(CLEAN),,$@) \
         -C $(subst $(CLEAN),,$@) clean

cleanall: clean $(EXAMPLECLEAN) 
