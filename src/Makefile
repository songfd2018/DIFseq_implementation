CC=g++

all: DIFseq_MCESM

# link
DIFseq_MCESM: DIFseq_MCESM.o rngstream.o
	$(CC) -fopenmp DIFseq_MCESM.o rngstream.o -o DIFseq_MCESM -L/Your/GSL/DIRECTORY/lib/ -lgsl -lgslcblas -lm
 
# compile
DIFseq_MCESM.o: DIFseq_MCESM.cpp
	$(CC) -c -fopenmp DIFseq_MCESM.cpp -w -std=c++11 -I/Your/GSL/DIRECTORY/include/

rngstream.o: rngstream.cpp
	$(CC) -c -fopenmp -lgomp rngstream.cpp

# clean all files
clean:
	rm -f *.o;
