# gnu compilers

# LAPACK and BLAS libraries
BLASLIB   = -L$(HOME)/LIB/LAPACK/3.8.0/BLAS -lblas_gfortran
LAPACKLIB = -L$(HOME)/LIB/LAPACK/3.8.0 -llapack_gfortran

# SuperLU library
SUPERLUDIR = $(HOME)/LIB/SuperLU/5.2.1
SUPERLULIB = -L$(SUPERLUDIR) -lsuperlu_gcc
SUPERLUINC = -I$(SUPERLUDIR)/SRC 

# ARPACK library
ARPACKDIR = $(HOME)/LIB/ARPACK
ARPACKLIB = -L$(ARPACKDIR) -larpack_gfortran
ARPACKINC = -I$(ARPACKDIR)/SRC

LIBS      = $(SUPERLULIB) $(ARPACKLIB) $(LAPACKLIB) $(BLASLIB) 
INCS      = $(SUPERLUINC) $(ARPACKINC)

CC        = gcc
CFLAGS    = -Wall -O2 -g

FC        = gfortran
FFLAGS    = -Wall -O2 -g

LOADER    = $(FC)
LDFLAGS   =

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

RM        = rm -f
