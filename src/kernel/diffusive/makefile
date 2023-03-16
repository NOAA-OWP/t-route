# compiler
FC := ${F90}

# compile flags
F90FLAGSGFORTRAN = -g -c -O2 -fPIC -fbounds-check
F90FLAGSINTEL = -free -g -c -O2 -fPIC -I${NETCDF}/include
ifeq ($(notdir ${FC}), ifort)
FCFLAGS = $(F90FLAGSINTEL)
else ifeq ($(notdir ${FC}), mpiifort)
FCFLAGS = $(F90FLAGSINTEL)
else
FCFLAGS = $(F90FLAGSGFORTRAN)
endif

diffusive.o: diffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pydiffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../src/troute-routing/troute/routing/fast_reach

clean:
	rm -f *.o *.mod
