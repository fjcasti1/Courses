SOURCES = quartile.cc 
OBJECTS = quartile.o
CC = icpc 
#CC = g++ 
OPT = 
LIBS = 
PARALLEL =
DEBUG = #-g
OTHER = -c -Wall
FFLAGS = $(OPT) $(OTHER) $(PARALLEL)

%.o: %.cc
	$(CC) $(FFLAGS) $<
%.o: %.c
	$(CC) $(FFLAGS) $<

all:  quartile

quartile: $(OBJECTS)
	$(CC) $(PARALLEL) $(OBJECTS) $(LIBS) -o quartile

clean:
	rm -rf $(OBJECTS) $(MODULES) 
	rm quartile
