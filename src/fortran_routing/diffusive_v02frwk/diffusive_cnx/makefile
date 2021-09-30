# compiler
FC := gfortran

# compile flags
FCFLAGS = -fbounds-check -g -c -O2 -fPIC

diffusive_cnx.o: diffusive_cnx.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive_cnx.o: pydiffusive_cnx.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../../src/python_routing_v02/troute/routing/fast_reach
	cp *.mod ../../../../src/python_routing_v02/troute/routing/fast_reach   
clean:
	rm -f *.o *.mod
