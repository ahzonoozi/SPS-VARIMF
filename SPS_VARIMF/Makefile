# Start of the Makefile
# Defining variables
#objects = main.o sps_setup.o 
#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
F90 = gfortran
#---------------------
# Intel
#F90 = ifort
#---------------------
#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
# GNU (version < 4.3.0)
F90FLAGS = -O -cpp
# GNU (version >= 4.3.0)
#F90FLAGS = -O3 -march=native -cpp

# Makefile
PROGS = main

COMMON = params.o sub_func.o  csp_out.o csp_gen.o Read_input.o\
         locate.o ssp_gen.o gal_mags.o imf_norm.o IGIMF_NORM.o\
         max_mecl.o max_mstar.o getspec.o filtr.o spec_wr.o spec_pagb.o\
         linterparr.o linterp.o sfh_gal.o spec_agb.o spec_wmb.o\
         imf_topheavy.o compute_indx.o

all : $(PROGS)

# Cleaning everything
clean :
	rm -rf *.o *.mod *.MOD *~ *.txt *.log fort.*

main : main.o $(COMMON)
	$(F90) -o main.exe main.o $(COMMON)

main.o : params.o sub_func.o

sps_setup.o :  params.o  sub_func.o

sub_func.o :  params.o

%.o : %.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<



