LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o) 
HEADERFILES =

all: $(EXAMPLENAME)

%.o: %.f90
	$(DIVER_FF) -c $(DIVER_FOPT) -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ)
	$(DIVER_FF) $(DIVER_FOPT) $(OBJ) $(LIBS) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
