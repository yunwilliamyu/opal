GDLDIR=ext/gdl-1.1/GDL

CC=gcc
CFLAGS= -O2 -Wall -Wshadow -I./ -I${GDLDIR}/include -c
OBJS = drawfrag.o fasta2vw.o fasta2skm.o fasta2sklearn.o

all: drawfrag fasta2vw fasta2skm fasta2sklearn

.c.o :
	$(CC) $(CFLAGS) $<

drawfrag : drawfrag.o
	$(CC) -o $@ drawfrag.o -lz -L${GDLDIR}/lib -lgdl -lm

fasta2vw : fasta2vw.o
	$(CC) -o $@ fasta2vw.o -lz -L${GDLDIR}/lib -lgdl -lm


fasta2skm : fasta2skm.o
	$(CC) -o $@ fasta2skm.o -lz -L${GDLDIR}/lib -lgdl -lm

fasta2sklearn: fasta2sklearn.o
	$(CC) -o $@ fasta2sklearn.o -lz -L${GDLDIR}/lib -lgdl -lm

clean :
	rm $(OBJS) drawfrag fasta2vw fasta2skm fasta2sklearn
