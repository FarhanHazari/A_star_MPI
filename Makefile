CC = mpicc
CFLAGS = -O3 -Wall -g -std=c11 -Wno-unused-function -Wno-deprecated-declarations
LDLIBS = -lm

a_star: a_star.o tools.o heap.o

.PHONY: clean
clean:
	rm -f *.o
	rm -f a_star
	rm -fr *.dSYM/
