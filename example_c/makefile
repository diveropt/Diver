LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o)
HEADERFILES = $(INC)/diver.h

all: $(EXAMPLENAME)

%.o: %.c
	$(DIVER_CC) -c $(DIVER_COPT) -std=c99 -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ) $(HEADERS)
	$(DIVER_FF) $(OBJ) $(LIBS) $(DIVER_MIXOPT_C) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
