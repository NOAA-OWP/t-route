FC := $(F90)

# compile flags
FFLAGS += -g -c -O2 -fPIC

diffusive.o: diffusive.f90
	$(FC) $(FFLAGS) -o $@ $<

pydiffusive.o: pydiffusive.f90
	$(FC) $(FFLAGS) -o $@ $<

install:
	cp *.o ../../../src/python_routing_v02/fast_reach

clean:
	rm -f *.o *.mod
