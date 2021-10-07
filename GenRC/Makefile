SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
CERNLIBS    := -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib

CXX       := g++
#CXXFLAGS  += -std=c++11 -O2 -Wall -fPIC $(ROOTCFLAGS)
CXXFLAGS  += -std=c++0x -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS) -lgfortran -Wl,-rpath,$(HAPRAD_CPP)/slib

INCLUDES  := -I$(HAPRAD_CPP)
#LIBS      := -L/$(HAPRAD_CPP)/slib -lTRadCor $(ROOTLIBS) -lMathMore -L/usr/lib64/cernlib/2006/lib/ $(CERNLIBS)
LIBS      := -L$(HAPRAD_CPP)/slib -lTRadCor $(ROOTLIBS) -lMathMore -L/opt/cern/pro/lib/ $(CERNLIBS)
##############################################################################
all: genRC

genRC:genRC.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

genRC.o:genRC.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf genRC genRC.o
