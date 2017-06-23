#! /bin/bash

CC=${1:-gcc}
CXX=${2:-g++}
CARCHFLAG=${3:--march=native}

CFLAGS="-O3 -Wall -Wno-unused-function -std=c99 -DBENCH $CARCHFLAG"
CXXFLAGS="-O3 -Wall -Wno-unused-function -std=c++11 -DBENCH $CARCHFLAG"

echo $CC$CARCHFLAG
make clean -s
make -s CFLAGS="$CFLAGS" CXXFLAGS="$CXXFLAGS"
for i in `seq 10`;
do
	echo $i;
	cat seq/pair.window.20.32.25k.txt | head -2000 | ./bench -i-;
done

