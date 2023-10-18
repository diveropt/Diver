LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o)
HEADERFILES = $(INC)/diver.hpp

all: $(EXAMPLENAME)

%.o: %.cpp
	$(DIVER_CC) -c $(DIVER_COPT) $(DIVER_CPPOPT) -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ) $(HEADERS)
	$(DIVER_FF) $(OBJ) $(LIBS) $(DIVER_MIXOPT_CPP) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
