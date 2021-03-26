# compiler
FC := gfortran

# compile flags
FCFLAGS = -g -c -O2 -fPIC

diffusive.o: main.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pymain.f90
	$(FC) $(FCFLAGS) -o $@ $<
    
install:
	cp *.o ../
    
clean:
	rm -f *.o *.mod
