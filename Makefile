
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

parasail: parasail/cpuid.c parasail/io.c parasail/matrix_lookup.c parasail/memory.c parasail/memory_avx2.c parasail/time.c parasail/sg_striped_avx2_256_16.c
	$(CC) $(CFLAGS) -c -o parasail/cpuid.o -I. parasail/cpuid.c
	$(CC) $(CFLAGS) -c -o parasail/io.o -I. parasail/io.c
	$(CC) $(CFLAGS) -c -o parasail/matrix_lookup.o -I. parasail/matrix_lookup.c
	$(CC) $(CFLAGS) -c -o parasail/memory.o -I. parasail/memory.c
	$(CC) $(CFLAGS) -c -o parasail/memory_avx2.o -I. parasail/memory_avx2.c
	$(CC) $(CFLAGS) -c -o parasail/time.o -I. parasail/time.c
	$(CC) $(CFLAGS) -c -o parasail/sg_striped_avx2_256_16.o -I. parasail/sg_striped_avx2_256_16.c


bench: bench.o aed.o alinear.o aaffine.o rlinear.o raffine.o gaba_wrap.o gaba_linear.o gaba_affine.o edlib.o ksw.o DB.o QV.o align.o seqan_wrap.o blast.o parasail/cpuid.o parasail/io.o parasail/matrix_lookup.o parasail/memory.o parasail/memory_avx2.o parasail/time.o parasail/sg_striped_avx2_256_16.o
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	rm -rf bench *.o
