# compiler
FC := gfortran

# compile flags
#FCFLAGS = -c -fdefault-real-8 -fno-align-commons -fbounds-check --free-form
FCFLAGS = -c -O2 -fPIC
# link flags
FLFLAGS = -static-gfortran -static-libgcc -no-defaultlibs -lgfortran -lgcc

.PHONY: reach reservoir all install

reach: mc_single_seg.o pymc_single_seg.o

mc_single_seg.o: MCsingleSegStime_f2py_NOLOOP.f90 precis.mod
	$(FC) $(FCFLAGS) -o $@ $<

pymc_single_seg.o: pyMCsingleSegStime_NoLoop.f90
	$(FC) $(FCFLAGS) -o $@ $<

reservoir: module_levelpool.o pymodule_levelpool.o

precis.mod: varPrecision.f90
	$(FC) $(FCFLAGS) $<

install:
	cp *.o ../../../../src/python_routing_v02/fast_reach

all: reach reservoir

%.o: %.F
	$(FC) $(FCFLAGS) -o $@ $<

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
