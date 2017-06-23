#! /bin/sh

for s in clang,clang++,-msse4.2 clang,clang++,-mavx clang,clang++,-march=native gcc,g++,-msse4.2 gcc,g++,-mavx gcc,g++,-march=native;
do
	./bench.sh `echo $s | sed 's/,/ /g'` > log.`echo $s | sed 's/,.*,//g'`.txt;
done

