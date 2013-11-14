# Makefile template for Diver 
# Author: Pat Scott, patscott@physics.mcgill.ca

# Until we sort out some configure magic, change 
# the options below to suit your system

# Define fortran compiler and options
FF=gfortran
CC=gcc
FOPT=-O3 -fPIC -ffree-line-length-none -ffixed-line-length-none -cpp
COPT=-O3 -fPIC
CPPOPT=-lstdc++
MIXOPT_C=
MIXOPT_CPP= -lstdc++
MODULE=-J

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

export FF FOPT CC COPT CPPOPT MIXOPT_C MIXOPT_CPP LIBNAME PREFIX LIB INC

SOURCEFILES = detypes deutils mutation crossover selection init converge posterior evidence io de cwrapper
OBJ = $(SOURCEFILES:%=$(BUILD)/%.o)

all: libdiver.a $(EXAMPLENAMES)

libdiver.a: makefile $(OBJ)
	$(AR) $(LIB)/$@ $(OBJ)

libdiver.so: makefile $(OBJ) 
	$(SHAR) -o $(LIB)/$@ $(OBJ)  
 
$(BUILD)/%.o: $(SOURCE)/%.f90
	cd $(BUILD); \
	$(FF) -c $(FOPT) -I$(INC) $(MODULE) $(INC) $<

$(BUILD)/%.o: $(SOURCE)/%.F90
	cd $(BUILD); \
	$(FF) -c $(FOPT) -I$(INC) $(MODULE) $(INC) $<

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
