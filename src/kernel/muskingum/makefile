# compiler
FC := $(F90)

# compile flags
#FCFLAGS = -c -fdefault-real-8 -fno-align-commons -fbounds-check --free-form
F90FLAGSGFORTRAN = -g -c -O2 -fPIC -lgfortran -lgcc -static-libgfortran -static-libgcc -nodefaultlibs
F90FLAGSINTEL = -free -w -c -O2 -fPIC -fconvert=big-endian -frecord-marker=4
ifeq ($(notdir ${FC}), ifort)
F90FLAGS = $(F90FLAGSINTEL)
else ifeq ($(notdir ${FC}), mpiifort)
F90FLAGS = $(F90FLAGSINTEL)
else
F90FLAGS = $(F90FLAGSGFORTRAN)
endif


.PHONY: reach all install
all: reach 

reach: mc_single_seg.o pymc_single_seg.o

mc_single_seg.o: MCsingleSegStime_f2py_NOLOOP.f90 precis.mod
	$(FC) $(F90FLAGS) -o $@ $<

pymc_single_seg.o: pyMCsingleSegStime_NoLoop.f90
	$(FC) $(F90FLAGS) -o $@ $<

precis.mod: varPrecision.f90
	$(FC) $(F90FLAGS) $<

install:
	cp *.o ../../../src/troute-routing/troute/routing/fast_reach

%.o: %.F
	$(FC) $(F90FLAGS) -o $@ $<

%.o: %.f90
	$(FC) $(F90FLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
