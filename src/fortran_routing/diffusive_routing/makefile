# compiler
FC := gfortran

# compile flags
FCFLAGS = -g -c -O2 -fPIC -fbounds-check

diffusive.o: diffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pydiffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../src/python_routing_v02/troute/routing/fast_reach

clean:
	rm -f *.o *.mod
