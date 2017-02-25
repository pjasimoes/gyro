# makefile for gyrosync in C++ program (PJSimoes, 2011)

CC=g++

CFLAGS=-Wall -c -O3
COPENMP=-fopenmp
CDEBUG=-g -Wall -c
CLIBS=
EXECUTABLE=-o gyro
SOURCES=gyro.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all:	
	$(CC) $(CFLAGS) $(COPENMP) $(SOURCES)
	$(CC) $(CLIBS) $(COPENMP) $(OBJECTS) $(EXECUTABLE)

single:	
	$(CC) $(CFLAGS) $(SOURCES)
	$(CC) $(CLIBS) $(OBJECTS) $(EXECUTABLE)

debug: 
	$(CC) $(CDEBUG) $(SOURCES)
	$(CC) $(CLIBS) $(OBJECTS) $(EXECUTABLE)
clean:
		rm -rf *.o

#
#g++ -Wall -I/usr/local/include drive.cpp -c	
#g++ drive.o -lgsl -lgslcblas -lm -o drive
