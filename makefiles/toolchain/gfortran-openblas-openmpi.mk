# Compiler and linker ==========================================================
FC                := gfortran
MPI_FC            := mpif90
FC                := $(MPI_FC)


# Fortran compiler flags =======================================================
FFLAGS            := -ffree-form -fno-range-check -ffree-line-length-none


# OpenMP option ================================================================
OPENMP            := -fopenmp


# Debug options ================================================================
DEBUG             := -O0 -g -fbacktrace -fcheck=bounds -Wall
DEBUG_FULL        := -Og -g3 -fbacktrace -fcheck=all -Wall -Wextra -Wunused-parameter -Wcharacter-truncation -Wsurprising -Waliasing -fimplicit-none -finit-real=nan -finit-integer=-999


# Fortran preprocessor =========================================================
FPP               := -cpp -DOPENBLAS


# .mod files related flags =====================================================
MODFLAG           := -J
IMOD              := -I


# Linking libraries provided by compiler environment ===========================
LDFLAGS           := -Wl,-rpath,$(OPENBLAS_ROOT)/lib
LIBS              := -L$(OPENBLAS_ROOT)/lib -lopenblas