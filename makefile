# makefile for gyrosync in C++ program (PJSimoes, 2011)

CC=g++

CFLAGS=-Wall -c -O3
COPENMP=-fopenmp
CDEBUG=-g -Wall -c
CLIBS=
EXECUTABLE=-o gyro
SOURCES=gyro.cpp
OBJECTS=$(SOURCES:.cpp=.o)

OPENMP_SUPPORTED := $(shell printf 'int main(){return 0;}\n' | $(CC) -x c++ -fopenmp - -o /tmp/gyro_openmp_test >/dev/null 2>&1 && echo yes || echo no)

ifeq ($(OPENMP_SUPPORTED),yes)
ALL_CFLAGS=$(CFLAGS) $(COPENMP)
ALL_LDFLAGS=$(CLIBS) $(COPENMP)
else
ALL_CFLAGS=$(CFLAGS)
ALL_LDFLAGS=$(CLIBS)
endif

all:	
	$(CC) $(ALL_CFLAGS) $(SOURCES)
	$(CC) $(ALL_LDFLAGS) $(OBJECTS) $(EXECUTABLE)

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
