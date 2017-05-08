CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c99 -DBENCH

CXX = g++
CXXFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c++11 -DBENCH

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

bench: bench.o aed.o alinear.o aaffine.o rlinear.o raffine.o gaba_wrap.o gaba_linear.o gaba_affine.o edlib.o
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	rm -rf bench *.o
