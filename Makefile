FC = mpiifort
# Fortran compiler (FC) (ifort: Intel compiler - recommended)

SRC_DIR = src
# Source code directory
OBJ_DIR = obj
# Objective file directory
MOD_DIR = mod
# Module file directory
BIN_DIR = bin
# Bin file directory

FFLG = -r8 -O3 -qopenmp -qmkl=parallel
# Fortran compiler flag
FINC = -I${MKLROOT}/include/fftw 
# Fortran include files
FLIB =  -L${MKLROOT}/lib/intel64 -liomp5 -lpthread \
		-lm -ldl
# Fortran libraries

SRC_F = $(wildcard $(SRC_DIR)/*.f90)
# List of source files.

# OBJ_C = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_C))
OBJ_F = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRC_F)) $(wildcard $(OBJ_DIR)/*.o)
# List of object files

ifneq ($(OBJ_DIR),)
  $(shell test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR))
endif
# If the obj folder not exists, make this directory

ifneq ($(MOD_DIR),)
  $(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
  FFLG+= -module $(MOD_DIR)
endif
# If the obj folder not exists, make this directory

ifneq ($(BIN_DIR),)
  $(shell test -d $(BIN_DIR) || mkdir -p $(BIN_DIR))
endif
# If the obj folder not exists, make this directory

EXE_F = qvortexgen \
	qvortexrun \
	dustyqvortgen \
	dustyqvortrun \
	test \
# Program scripts in the f90 folder (ADD/REMOVE THE PROGRAM LISTS HERE)
# Using 'make [program_name]' will create the executable program file in the bin folder
# 'make new' will wipe out all obj, mod and exec files and re-compile the first exe_f file

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(FC) -c $< -o $@  $(FFLG) $(FINC) $(FLIB)
# Create an object file from c source files
# Module files go to the MOD_DIR directory

all: $(OBJ_F)

.PHONY: clean new

$(EXE_F): $(OBJ_F)
	$(FC) ./f90/$@.f90 $^ -o $(BIN_DIR)/$@_exec $(FFLG) $(FINC) $(FLIB)
# Produce an executable file labeled by _exec

clean:
	rm -f ./*.dat
	find ./$(OBJ_DIR)/ ! -name 'fm*' -exec rm -f {} \;
	find ./$(MOD_DIR)/ ! -name 'fm*' -exec rm -f {} \;
	rm -f ./$(BIN_DIR)/*
# Remove all objects, modules and executable files
# Due to compilation time, the fm objects and modules are not deleted
# If you want to remove them, manually delete them

new:
	make clean
	make
# Clean and produce an executable file from the first program script in EXE_F

include Makefile.dependencies
# Reference to module dependencies
