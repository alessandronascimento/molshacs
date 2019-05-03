#!/bin/bash

CC=g++
CFLAGS=-c -O3 -ffast-math -DNLOPT -w
LDFLAGS= -lgsl -lgslcblas -lm -lz -lnlopt_cxx -static
SOURCES=CORREL.cpp  Gaussian.cpp Grid.cpp  main.cpp  Minimizer2.cpp  Mol2.cpp  Parser.cpp  RunEngine.cpp  Writer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MolShaCS.x
EXE_PATH=/home/asn/bin

all: $(SOURCES) $(EXECUTABLE) install 
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) $(EXECUTABLE)

install:
	cp $(EXECUTABLE) $(EXE_PATH)
