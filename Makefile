
FLAGS = -O3 -Wall -Wno-unused-function -march=native -DBENCH -DCOUNT_CELLS -I.

CC = gcc
CFLAGS = $(FLAGS) -std=c99

CXX = g++
CXXFLAGS = $(FLAGS) -std=c++14

all: bench

.c.o:
	$(CC) -c $(CFLAGS) $<

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $<

gaba_linear.o: gaba/gaba.c
	$(CC) -c -o $@ $(CFLAGS) -DMODEL=LINEAR -DSUFFIX $<

gaba_affine.o: gaba/gaba.c
	$(CC) -c -o $@ $(CFLAGS) -DMODEL=AFFINE -DSUFFIX $<

gaba_wrap.o: gaba/gaba_wrap.c
	$(CC) -c $(CFLAGS) $<

seqan_wrap.o: seqan_wrap.cc
	$(CXX) -c $(CXXFLAGS) $<

bench: bench.o aed.o alinear.o aaffine.o rlinear.o raffine.o gaba_wrap.o gaba_linear.o gaba_affine.o edlib.o ksw.o DB.o QV.o align.o seqan_wrap.o blast.o
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	rm -rf bench *.o
