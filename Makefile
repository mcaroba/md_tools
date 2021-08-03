# Copyright (c) 2021 by Miguel A. Caro

SHELL = /bin/sh

# User-modifiable variables
F90=gfortran

# Default locations for various files
BUILD_DIR=build
BIN_DIR=bin
INC_DIR=include
LIB_DIR=lib

# Edit F90_OPTS, F90_MOD_DIR_OPT and LIBS as required
#########################################################
ifeq ($(F90),gfortran)
F90_OPTS=-fPIC -O3 -fopenmp
#
# For debugging:
#F90_OPTS=-fPIC -O3 -fcheck=bounds -g -fcheck=all -Wall
#
# Preprocessing is not currently used, but will be needed
# if I add MPI support.
PP=-cpp
F90_MOD_DIR_OPT=-J
LIBS=
endif

# Do not change anything below this line
##########################################################

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := get_steinhardt_ortho slice_xyz

SRC := potentials.f90 neighbors.f90 analyze.f90
OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))
PROG := $(addprefix $(BIN_DIR)/,$(PROGRAMS))

.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: default all programs clean deepclean libmdtools

default: libmdtools programs

all: default

clean:
	rm -rf $(OBJ) $(INC_DIR)/*.mod $(PROG) 

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)

libmdtools: $(OBJ) ${LIB_DIR}
	ar scr $(LIB_DIR)/libmdtools.a $(OBJ)

$(BIN_DIR)/%: src/%.f90 $(OBJ) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o: src/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR): ${INC_DIR}
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(INC_DIR):
	mkdir -p $@

$(LIB_DIR):
	mkdir -p $@
