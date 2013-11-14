LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o)
HEADERFILES = $(INC)/diver.hpp

all: $(EXAMPLENAME)

%.o: %.cpp
	$(CC) -c $(COPT) $(CPPOPT) -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ) $(HEADERS)
	$(FF) $(MIXOPT_CPP) $(OBJ) $(LIBS) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
