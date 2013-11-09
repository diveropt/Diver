# Makefile template for DEvoPack/Fortran 
# Author: Pat Scott, patscott@physics.mcgill.ca

# Until we sort out some configure magic, change 
# the options below to suit your system

# Define fortran compiler and options
FF=gfortran
CC=gcc
FOPT=-O3 -fPIC -ffree-line-length-none -ffixed-line-length-none
COPT=-O3 -fPIC
CPPOPT=-lstdc++
MIXOPT_C=
MIXOPT_CPP= -lstdc++
MODULE=-J

CLEAN = _clean
AR = ar r  
SHAR = ld -shared

LIBNAME = devo
EXAMPLENAMES = example_f example_c example_cpp
EXAMPLECLEAN = $(EXAMPLENAMES:%=%$(CLEAN))
PREFIX = devo

DEVO_ROOT = $(PWD)
SOURCE = ${DEVO_ROOT}/source
INC = ${DEVO_ROOT}/include
BUILD = ${DEVO_ROOT}/build
LIB = ${DEVO_ROOT}/lib

export FF FOPT CC COPT CPPOPT MIXOPT_C MIXOPT_CPP LIBNAME PREFIX LIB INC

SOURCEFILES = detypes deutils mutation crossover selection init converge posterior evidence io de cwrapper
OBJ = $(SOURCEFILES:%=$(BUILD)/%.o)

all: libdevo.a $(EXAMPLENAMES)

libdevo.a: makefile $(OBJ)
	$(AR) $(LIB)/$@ $(OBJ)

libdevo.so: makefile $(OBJ) 
	$(SHAR) -o $(LIB)/$@ $(OBJ)  
 
$(BUILD)/%.o: $(SOURCE)/%.f90
	cd $(BUILD); \
	$(FF) -c $(FOPT) -I$(INC) $(MODULE) $(INC) $<

$(BUILD)/%.o: $(SOURCE)/%.F90
	cd $(BUILD); \
	$(FF) -c $(FOPT) -I$(INC) $(MODULE) $(INC) $<

$(EXAMPLENAMES): libdevo.a
	make EXAMPLENAME=$@ -C $@

clean:
	rm -f $(LIB)/*.a $(LIB)/*.so; \
	cd $(BUILD); rm -f *.o *.mod; \
	cd $(INC); rm -f *.mod

$(EXAMPLECLEAN):
	make EXAMPLENAME=$(subst $(CLEAN),,$@) \
         -C $(subst $(CLEAN),,$@) clean

cleanall: clean $(EXAMPLECLEAN) 
