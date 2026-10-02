# Validation harness for HAPRAD_CPP.
#
#   make            build rc_point; build haprad2_point too if CERN_LIB is set
#   make compare    build, then run the default grid through both codes
#   make clean
#
# rc_point       -- the C++ implementation (this repository)
# haprad2_point  -- the original HAPRAD 2.0 Fortran in ../haprad2, which is the
#                   code the C++ was ported from. It needs PDFLIB (PDFSET,
#                   STRUCTM) for its PDF x FF semi-inclusive model, so it is
#                   built only when CERN_LIB points at a libpdflib804.

.PHONY: all clean compare
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

TOP        := ..
HAPRAD2    := $(TOP)/haprad2
BIN        := bin
OBJ        := .obj

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS:= $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS   := $(shell $(ROOTCONFIG) --libs)

CXXFLAGS   := -g -O2 -Wall -I$(TOP) $(ROOTCFLAGS)
LDFLAGS    := -g -O2 $(ROOTLDFLAGS) -Wl,-rpath,$(abspath $(TOP))/slib
LIBS       := -L$(TOP)/slib -lTRadCor $(ROOTLIBS) -lMathMore

FC         := gfortran
FFLAGS     := -std=legacy -fno-automatic -ffixed-line-length-none -fdollar-ok \
              -fno-second-underscore -w -I$(HAPRAD2)/include

# Everything from haprad2 except its own main programs and parameter setup.
F2SRC      := ihaprad.f haprad_utils.f exclusive_model.f h3.f h4.f \
              pkhff.f init_pdf.f semi_inclusive_model.f
F2OBJ      := $(addprefix $(OBJ)/,$(F2SRC:.f=.o))

TARGETS    := $(BIN)/rc_point $(BIN)/phi_convention $(BIN)/make_toy_table $(BIN)/make_model_table $(BIN)/make_toy_events \
              $(BIN)/make_weighted_events
ifdef CERN_LIB
TARGETS    += $(BIN)/haprad2_point $(BIN)/born_harmonics
endif

# For tabulating HAPRAD 2.0's model (PLAN.md 2.3): the same objects, but with
# the unconditional `stop` at semi_inclusive_model.f:178 removed, so the cap
# the authors wrote immediately after it takes effect instead of aborting the
# run. The reference used by haprad2_point is not touched.
F2TABOBJ   := $(filter-out $(OBJ)/semi_inclusive_model.o,$(F2OBJ)) $(OBJ)/semi_inclusive_model_nostop.o

# Both codes open their data files from the current working directory:
# pi_n_maid.dat (MAID exclusive model) and the *.grid PKH fragmentation
# functions. Link them in rather than copying 2.4 MB around.
DATA := pi_n_maid.dat plo.grid pnlo.grid klo.grid knlo.grid hlo.grid hnlo.grid
DATALINKS := $(addprefix data/,$(DATA))

data:
	@mkdir -p $@

$(DATALINKS): | data
	@ln -sf ../$(HAPRAD2)/$(@F) $@

all: $(TARGETS) $(DATALINKS)
ifndef CERN_LIB
	@echo ""
	@echo "NOTE: CERN_LIB is not set, so the haprad2 reference was NOT built."
	@echo "      Only rc_point (self-consistency checks) is available."
	@echo "      Build PDFLIB, then:  export CERN_LIB=\$$HOME/cernlib/lib && make"
	@echo ""
endif

$(BIN) $(OBJ):
	@mkdir -p $@

$(BIN)/rc_point: src/rc_point.cxx | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@

$(BIN)/make_model_table: src/make_model_table.cxx | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@

$(BIN)/make_toy_events: src/make_toy_events.cxx | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(ROOTLIBS) -o $@

$(BIN)/make_weighted_events: src/make_weighted_events.cxx | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(ROOTLIBS) -o $@

$(BIN)/make_toy_table: src/make_toy_table.cxx | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(ROOTLIBS) -o $@

$(BIN)/phi_convention: src/phi_convention.cxx | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@

$(OBJ)/%.o: $(HAPRAD2)/%.f | $(OBJ)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ)/semi_inclusive_model_nostop.o: $(HAPRAD2)/semi_inclusive_model.f | $(OBJ)
	sed 's/^      if (abs(4.d0 \* m_cos_phi) .gt. 0.9d0) stop$$/C     stop removed for tabulation; the cap below applies/' \
	    $< > $(OBJ)/semi_inclusive_model_nostop.f
	@if grep -qE '^ +if .* stop *$$' $(OBJ)/semi_inclusive_model_nostop.f; then \
	    echo "ERROR: the stop in semi_inclusive_model.f was not removed"; exit 1; fi
	$(FC) $(FFLAGS) -c $(OBJ)/semi_inclusive_model_nostop.f -o $@

$(BIN)/born_harmonics: src/born_harmonics.f $(F2TABOBJ) | $(BIN)
	$(FC) $(FFLAGS) $^ -L$(CERN_LIB) -Wl,-rpath,$(CERN_LIB) -lpdflib804 -o $@

$(BIN)/haprad2_point: src/haprad2_point.f $(F2OBJ) | $(BIN)
	$(FC) $(FFLAGS) $^ -L$(CERN_LIB) -Wl,-rpath,$(CERN_LIB) -lpdflib804 -o $@

compare: all
	./compare.py

clean:
	rm -rf $(BIN) $(OBJ)
