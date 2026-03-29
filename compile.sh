#!/bin/bash
#
# Compilation script for CRS model using gfortran
#

# Compiler and flags
FC=gfortran
FFLAGS="-O2  -ffree-form -ffree-line-length-none -ffpe-summary=none"
DEBUGFLAGS="-g -fcheck=all -fbacktrace -Wall -Wextra"

# Output executable name
EXE=crs_model

# Check for debug mode
if [ "$1" == "debug" ]; then
    FFLAGS="$DEBUGFLAGS -fopenmp"
    EXE=crs_model_debug
    echo "Compiling in DEBUG mode..."
else
    echo "Compiling in RELEASE mode..."
fi

# Clean old files
rm -f *.o *.mod 
rm -f $EXE

# Source files in dependency order
# 1. Modules first (crs_mod contains datatype, sizes, globals, pickmoments, myseed)
# 2. Then subroutines that use those modules
# 3. Finally the main program

echo "Compiling modules..."
$FC $FFLAGS -c crs_mod.f90

echo "Compiling model subroutines..."
$FC $FFLAGS -c modelfunctions.f90
$FC $FFLAGS -c shock_draw.f90
$FC $FFLAGS -c simmodel.f90
$FC $FFLAGS -c makemoments.f90
$FC $FFLAGS -c crs_solve.f90

echo "Compiling main program..."
$FC $FFLAGS -c evaluate_crs.f90

echo "Linking..."
$FC $FFLAGS -o $EXE \
    crs_mod.o \
    modelfunctions.o \
    shock_draw.o \
    simmodel.o \
    makemoments.o \
    crs_solve.o \
    evaluate_crs.o

if [ $? -eq 0 ]; then
    echo "Compilation successful! Executable: $EXE"
    echo ""
    echo "Run with: ./$EXE"
    echo "Set threads with: export OMP_NUM_THREADS=N"
else
    echo "Compilation failed!"
    exit 1
fi
