
# compile flags
FCFLAGS = -g -c -O2 -fPIC

diffusive.o: diffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pydiffusive.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../src/python_routing_v02/fast_reach

clean:
	rm -f *.o *.mod
