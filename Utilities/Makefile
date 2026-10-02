#SHELL = /bin/bash

.DELETE_ON_ERROR:

MKDIR_P := mkdir -p

BINDIR := ./bin
SRCDIR := ./src

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS := $(shell $(ROOTCONFIG) --libs)
# CERNLIB is not required: libTRadCor.so pulls in no CERNLIB symbols unless it
# was built with USE_CERNLIB=1. Link against it only when asked for, so that a
# machine without CERNLIB (and CERN_LIB unset) still builds these programs.
CERNLIBS :=
CERNLIB_LIBS ?= -lpdflib804
ifeq ($(USE_CERNLIB),1)
    CERNLIBS := -L$(CERN_LIB) $(CERNLIB_LIBS)
endif

CXX := c++ # -D DEBUG
LD := c++

CXXFLAGS := -g -O2 -Wall -fPIC $(ROOTCFLAGS)
LDFLAGS := -g -O2 $(ROOTLDFLAGS) -lgfortran -Wl,-rpath,$(HAPRAD_CPP)/slib

INCLUDES := -I./include -I$(HAPRAD_CPP)
LIBS := -L$(HAPRAD_CPP)/slib -lTRadCor $(ROOTLIBS) -lMathMore $(CERNLIBS)

PROG := GetBinning GetCentroids FitPhiPQ GetRC MakePhiTable MakeRCGrid ApplyRC
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
