#!/bin/sh

CC=gcc
OBJECTS = mersenne.o 
SRCDIR = ./src-clstr

algorithm: $(OBJECTS)
	$(CC) $(OBJECTS) algorithm.c -Wall -Wextra -Wformat -Wformat-security -o algorithm -lpthread -lm -O3 -pg 

mersenne.o: mersenne.c
	$(CC) mersenne.c -c -lm -Wall -Wextra -O3

clean:
	rm algorithm  mersenne.o

