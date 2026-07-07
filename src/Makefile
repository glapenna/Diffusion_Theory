
BASEDIR=${HOME}/diffusion
EXEDIR=$(BASEDIR)/bin
#EXENAME=mdan.npbc
EXENAME=mdan_gpu.npbc
CC=pgcc
#CC=gcc
FC=pgfortran
#FC=gfortran

# standard
#CFLAGS= -DHYDRO_OSEEN -O2 -acc=gpu -Minfo=all # -gpu=cc80
#CFLAGS= -DHYDRO_OSEEN -g -acc=gpu -gpu=cc80
CFLAGS= -DHYDRO_OSEEN -O2 -acc=gpu -gpu=cc80
#CFLAGS= -DHYDRO_OSEEN -O2
#CFLAGS= -DHYDRO_OSEEN -ggdb -fsanitize=address -fno-omit-frame-pointer
#LIBFLAGS = -lm -O2 -acc=gpu -Minfo=all # -gpu=cc80 -lc
#LIBFLAGS = -lm -O2 -acc=gpu -gpu=cc80
#LIBFLAGS = -lm -O2
LIBFLAGS = -lm

OBJ= main.o stream.o doit.o rdprm.o rdprm_his2.o rdprm_his3.o rdprm_pdb.o \
        coord.o modsel.o pdbgen.o util.o \
        bond.o angle.o dihe.o gr.o geom.o sputinbox.o \
        utilrot.o inertdir.o nrutil.o stats.o \
        orderb.o geomq.o density1.o strlib.o gyr.o \
	table.o bead.o mcarig_m1.o acff.o \
        rotdiff.o  mca_p2p1.o \
        complex.o svd.o fft.o group.o asa.o pucker.o \
	static.o rms.o orient.o weight.o dipol.o ortho.o pers.o map.o \
        asaf.o D2mn.o \
        hbond.o

install_double: mdanal
	mv mdanal $(EXEDIR)/$(EXENAME)

install_single: mdanal
	mv mdanal $(EXEDIR)/$(EXENAME)s

mdanal:	$(OBJ) 
	$(CC) $(CFLAGS) -o mdanal $(OBJ) $(LIBFFLAGS) $(LIBFLAGS)

#	$(FC) $(CFLAGS) -o mdanal $(OBJ) $(LIBFFLAGS) $(LIBFLAGS)

$(OBJ): head.h prm.h complex.h

mcarig_m1.o: head.h prm.h complex.h fbb.c fbb0.c hydrod.c

mca_p2p1.o: head.h prm.h complex.h fbb.c hydrod.c

rotdiff.o: head.h prm.h complex.h hydrod.c

D2mn.o: head.h nr.h complex.h

acff.o: head.h nr.h complex.h acff.c

clean:
	/bin/rm *.o
