#SHELL = /bin/bash

.DELETE_ON_ERROR:

MKDIR_P := mkdir -p

BINDIR := ./bin
SRCDIR := ./src

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS := $(shell $(ROOTCONFIG) --libs)
CERNLIBS := -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib

CXX := c++ # -D DEBUG
LD := c++

CXXFLAGS := -g -O2 -Wall -fPIC $(ROOTCFLAGS)
LDFLAGS := -g -O2 $(ROOTLDFLAGS) -lgfortran -Wl,-rpath,$(HAPRAD_CPP)/slib

INCLUDES := -I./include -I$(HAPRAD_CPP)
LIBS := -L$(HAPRAD_CPP)/slib -lTRadCor $(ROOTLIBS) -lMathMore -L$(CERN_LIB) $(CERNLIBS)

PROG := GetBinning GetCentroids FitPhiPQ GetRC
LIST := $(addprefix $(BINDIR)/, $(PROG))

.PHONY: all clean

all: $(LIST)

$(BINDIR)/%:$(BINDIR)/%.o
	@echo "Doing application" $@
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

$(BINDIR)/%.o:$(SRCDIR)/%.cxx
	$(MKDIR_P) $(BINDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf $(LIST:%=%.o) $(LIST)
