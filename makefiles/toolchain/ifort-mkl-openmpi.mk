# Compiler and linker ==========================================================
FC                := ifort
MPI_FC            := mpif90
FC                := $(MPI_FC)


# Fortran compiler flags =======================================================
FFLAGS            := 


# OpenMP option ================================================================
OPENMP            := -qopenmp


# Debug options ================================================================
DEBUG             := -O0 -traceback -check bounds
DEBUG_FULL        := -O0 -g -traceback -check bounds -check uninit -check pointers -ftrapuv -warn interfaces -gen-interfaces


# Fortran preprocessor =========================================================
FPP               := -fpp -DMKL


# .mod files related flags =====================================================
MODFLAG           := -module 
IMOD              := -I


# Linking libraries provided by compiler environment ===========================
LDFLAGS           := -Wl,-rpath,$(MKLROOT)/lib/intel64
LIBS              := -qmkl=parallel
