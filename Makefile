SHELL :=/bin/bash

#================================

CC=gcc-8
CFLAGS=-I. -O3 -fopenmp
DEPS = pagerank.h
OBJ = pagerank.o  

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pagerank: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) 

clean:
	$(RM) *.o *~
	
