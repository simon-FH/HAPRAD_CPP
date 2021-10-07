SHELL = /bin/bash

.DELETE_ON_ERROR:

.SECONDARY: $(DICT_CLASS)
.PHONY: all lib clean distclean checkdirs


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTLIBDIR  := $(shell $(ROOTCONFIG) --libdir)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)
ROOTCINT    := rootcint

ifeq "$(DEBUG)" ""
CXXFLAGS += -O2 -fPIC
else
CXXFLAGS += -g -fPIC
endif

#ifeq ($(findstring Linux,$(OS_NAME)),Linux)
#CXX       := g++  -fPIC 
#CXXFLAGS  += -Wall $(ROOTCFLAGS) $(SET_DEBUG) \
#    ./configure--enable-shared_dependencies
#endif

LD        = g++ -fPIC 
LDFLAGS   = -O2 $(ROOTLDFLAGS) -lgfortran 
SOFLAGS   = -Wl,-soname,$(notdir $@) -shared

INCLUDES  := -I$(ROOTINCDIR)
LOCALINCLUDES := -I.
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

FSRC_OBJ   := $(addprefix $(OBJ_DIR)/,$(patsubst %.f,%.o,$(wildcard *.f)))

SH_LIB     := libTRadCor.so

CERNLIBS   := -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib
FLIBS      := $(CERNLIBS)

FCC        := gfortran

F77OPT     := -fPIC -DLinux -ffixed-line-length-none -fdollar-ok \
              -fno-second-underscore

##############################################################################

$(OBJ_DIR)/%.o: %.f
	$(FCC) $(F77OPT) -c $< -o $@

##############################################################################

all: lib

lib: checkdirs $(SLIB_DIR)/$(SH_LIB)

include Makefile_depends


$(SLIB_DIR)/$(SH_LIB): $(SRC_OBJ) $(DICT_OBJ) $(FSRC_OBJ)
	$(LD) $(SOFLAGS) $(LDFLAGS)  $^ $(LIBS) -o $@


$(DICT_DIR)/%Dict.cxx: %.h
ifneq (,$(filter %.o,$(MAKECMDGOALS)))
	@test -d $(DICT_DIR) || mkdir -p $(DICT_DIR)
endif
	$(ROOTCINT) -f $@ -c $(<F) $(<F:%.h=%LinkDef.h)


$(OBJ_DIR)/%Dict.o: $(DICT_DIR)/%Dict.cxx
ifneq (,$(filter %.o,$(MAKECMDGOALS)))
	@test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)
endif
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LOCALINCLUDES) -c $< -o $@


$(OBJ_DIR)/%.o: %.cxx
ifneq (,$(filter %.o,$(MAKECMDGOALS)))
	@test -d $(DEP_DIR) || mkdir -p $(DEP_DIR)
	@test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)
endif
	@$(call make-depend,$<,$@,$(@F:.o=.d))
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LOCALINCLUDES) -c $< -o $@


checkdirs: $(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR)

$(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(OBJ_DIR) $(DICT_DIR)
	@rm -rf $(SLIB_DIR) $(DEP_DIR)
