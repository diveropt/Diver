LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o)
HEADERFILES = $(INC)/diver.h

all: $(EXAMPLENAME)

%.o: %.c
	$(DIVER_CC) -c $(DIVER_COPT) -std=c99 -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ) $(HEADERS)
	$(DIVER_FF) $(DIVER_MIXOPT_C) $(OBJ) $(LIBS) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
