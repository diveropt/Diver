LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o)
HEADERFILES = $(INC)/devo.h

all: $(EXAMPLENAME)

%.o: %.cpp
	$(CC) -c $(COPT) -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ) $(HEADERS)
	$(FF) $(MIXOPT) $(OBJ) $(LIBS) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
