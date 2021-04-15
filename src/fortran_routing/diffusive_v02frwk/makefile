# compiler
FC := gfortran

# compile flags
FCFLAGS = -g -c -O2 -fPIC

diffusive.o: maindiff.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pymaindiff.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../src/python_routing_v02/fast_reach

clean:
	rm -f *.o *.mod