SHELL = /bin/bash

# ifndef ROOTSYS
#     $(error "ERROR: please, set ROOT")
# endif

.DELETE_ON_ERROR:

.SECONDARY: $(DICT_CLASS)
.PHONY: all lib clean distclean checkdirs

ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cxx --cflags)
ROOTCFLAGS  := $(ROOTCFLAGS) -fPIC
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --cxx --ldflags)
ROOTLDFLAGS := $(ROOTLDFLAGS) -lgfortran
ROOTCINT    := rootcint

SOFLAGS = -Wl,-soname,$(notdir $@) -shared

LOCALINCLUDES := -I.

ROOTLIBDIR  := $(shell $(ROOTCONFIG) --libdir)
LIBS      := -L$(ROOTLIBDIR) -lMathMore

SLIB_DIR := slib
OBJ_DIR  := .obj
DICT_DIR := .dict
DEP_DIR  := .dep

SRC_CODE   := $(wildcard *.cxx)
SRC_CLASS  := TRadCor.cxx
SRC_DEP    := $(addprefix $(DEP_DIR)/,$(SRC_CODE:.cxx=.d))
SRC_OBJ    := $(addprefix $(OBJ_DIR)/,$(SRC_CODE:.cxx=.o))
DICT_CLASS := $(addprefix $(DICT_DIR)/,$(SRC_CLASS:.cxx=Dict.cxx))
DICT_OBJ   := $(addprefix $(OBJ_DIR)/,$(SRC_CLASS:.cxx=Dict.o))
FSRC_ALL   := $(wildcard *.f)

# init_pdf.f is nothing but a wrapper around the CERNLIB/PDFLIB entry points
# PDFSET and STRUCTM. No C++ source in this project calls it, so building it in
# only serves to drag two unresolvable symbols into libTRadCor.so, which then
# breaks the link of every program that uses the library. It is therefore left
# out by default. Build with `make USE_CERNLIB=1` (and CERN_LIB pointing at the
# CERNLIB libdir) to compile and link it -- needed only if the PDF x FF
# semi-inclusive model in TSemiInclusiveModel.cxx is ever restored.
#
# Only PDFLIB is actually required: this project ships its own dfint
# (THapradUtils.cxx) and integrates via GSL/ROOT, so packlib/kernlib/mathlib are
# not needed. Override CERNLIB_LIBS if your CERNLIB install requires more.
CERNLIB_LIBS ?= -lpdflib804
ifeq ($(USE_CERNLIB),1)
    LIBS := $(LIBS) -L$(CERN_LIB) -Wl,-rpath,$(CERN_LIB) $(CERNLIB_LIBS)
else
    FSRC_ALL := $(filter-out init_pdf.f,$(FSRC_ALL))
endif

FSRC_OBJ   := $(addprefix $(OBJ_DIR)/,$(patsubst %.f,%.o,$(FSRC_ALL)))

SH_LIB     := libTRadCor.so

FCC        := gfortran

F77OPT     := -std=legacy -fPIC -DLinux -ffixed-line-length-none -fdollar-ok -fno-second-underscore

$(OBJ_DIR)/%.o: %.f
	$(FCC) $(F77OPT) -c $< -o $@

all: lib

lib: checkdirs $(SLIB_DIR)/$(SH_LIB)

include Makefile_depends

$(SLIB_DIR)/$(SH_LIB): $(SRC_OBJ) $(DICT_OBJ) $(FSRC_OBJ)
	$(ROOTLDFLAGS) $(SOFLAGS) $^ $(LIBS) -o $@

$(DICT_DIR)/%Dict.cxx: %.h
	$(ROOTCINT) -f $@ -c $(<F) $(<F:%.h=%LinkDef.h)

$(OBJ_DIR)/%Dict.o: $(DICT_DIR)/%Dict.cxx
	$(ROOTCFLAGS) $(LOCALINCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: %.cxx
	@$(call make-depend,$<,$@,$(@F:.o=.d))
	$(ROOTCFLAGS) $(LOCALINCLUDES) -c $< -o $@

checkdirs: $(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR)

$(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR):
	@mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(DICT_DIR) $(SLIB_DIR) $(DEP_DIR)
