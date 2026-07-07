COMP= gcc
FLIB=
LFLAGS= 
#for GROMACS
#CFLAGS= -c -I/usr/share/gromacs-3.3/include
CFLAGS= -c -g -I/usr/include/tirpc
#LIBFLAGS= -L/usr/local/gromacs/lib -lgmx -lXt -lXm -lm
LIBFLAGS= -g -lm -ltirpc

EXEDIR=$(HOME)/bin
EXENAME=bin2xtc

OBJ= bin2xtc.o libxdr.o


install: all
	mv $(EXENAME) $(EXEDIR)

all:	$(OBJ)
	$(COMP) $(LFLAGS) -o $(EXENAME) $(OBJ) $(LIBFLAGS)

$(OBJ):

bin2xtc.o: bin2xtc.c
	$(COMP) $(CFLAGS) bin2xtc.c

libxdr.o: libxdr.c
	$(COMP) $(CFLAGS) libxdr.c

clean:
	/bin/rm *.o
