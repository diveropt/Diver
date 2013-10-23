LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o) 
HEADERFILES =

all: $(EXAMPLENAME)

%.o: %.f90
	$(FF) -c $(FOPT) -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ)
	$(FF) $(FOPT) $(OBJ) $(LIBS) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
