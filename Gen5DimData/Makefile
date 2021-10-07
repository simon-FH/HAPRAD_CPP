SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)

CXX       := g++
#CXXFLAGS  += -std=c++11 -O2 -Wall -fPIC $(ROOTCFLAGS)
CXXFLAGS  += -std=c++0x -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS)

INCLUDES  := -I/$(ROOTINCDIR)
LIBS      := $(ROOTLIBS)

##############################################################################
all: 5dim_acc_corr

5dim_acc_corr:5dim_acc_corr.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

5dim_acc_corr.o:5dim_acc_corr.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf 5dim_acc_corr 5dim_acc_corr.o
