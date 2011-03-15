#!/usr/bin/make
CC = mpicc
CFLAGS=-g -Wall
LDFLAGS=

all: life mkinit
mkinit: mkinit.c
	gcc $(CFLAGS) $? $(LDFLAGS) -o $@
clean:
	rm -f life mkinit *.o
