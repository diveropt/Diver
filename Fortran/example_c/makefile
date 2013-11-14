LIBS = -L$(LIB) -l$(LIBNAME)
NAME = ../$(PREFIX).$(EXAMPLENAME)
SOURCEFILES = $(EXAMPLENAME)
OBJ = $(SOURCEFILES:%=%.o)
HEADERFILES = $(INC)/diver.h

all: $(EXAMPLENAME)

%.o: %.c
	$(CC) -c $(COPT) -std=c99 -I$(INC) $<
 
$(EXAMPLENAME): $(OBJ) $(HEADERS)
	$(FF) $(MIXOPT_C) $(OBJ) $(LIBS) -o $(NAME)

clean:
	rm -f *.o *.mod $(NAME) 
