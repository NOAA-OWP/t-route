# compiler
FC := ${F90}

# compile flags
ifeq (${FC},gfortran)
FCFLAGS = -g -c -O2 -fPIC -fbounds-check
else
FCFLAGS = -free -g -c -O2 -fPIC -I${NETCDF}/include
endif

diffusive.o: diffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pydiffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../src/troute-routing/troute/routing/fast_reach

clean:
	rm -f *.o *.mod
