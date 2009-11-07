#!/usr/bin/make
CC = mpicc
CFLAGS=-g -Wall
LDFLAGS=

all: life mkinit
clean:
	rm -f life mkinit
