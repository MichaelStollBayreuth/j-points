all: j-points

j-points-fast: j-points-1.1.o j-sift-1.0-i386.s
	gcc j-points-1.1.o j-sift-1.0-i386.s -o j-points-fast \
	    -lgmp -lgcc -lc -lm

j-points: j-points-1.1.o j-sift-1.0.o
	gcc j-points-1.1.o j-sift-1.0.o -o j-points \
	    -lgmp -lgcc -lc -lm

j-points-1.1.o: j-points-1.1.c j-points.h
	gcc j-points-1.1.c -Wall -c -o j-points-1.1.o -O -fomit-frame-pointer

j-sift-1.0.o: j-sift-1.0.c j-points.h
	gcc j-sift-1.0.c -Wall -c -o j-sift-1.0.o -O -fomit-frame-pointer


