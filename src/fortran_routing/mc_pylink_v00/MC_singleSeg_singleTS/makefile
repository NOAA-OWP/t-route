# compiler
FC := $(F90)

# compile flags
#FFLAGS = -c -fdefault-real-8 -fno-align-commons -fbounds-check --free-form
FFLAGS += -g -c -O2 -fPIC
# link flags
# TODO: find out where these should be used (they are not, at the moment...)
LDFLAGS += -static-gfortran -static-libgcc -no-defaultlibs -lgfortran -lgcc
VPATH = ../Reservoir_singleTS

.PHONY: reach reservoir all install
all: reach reservoir

reach: mc_single_seg.o pymc_single_seg.o

mc_single_seg.o: MCsingleSegStime_f2py_NOLOOP.f90 precis.mod
	$(FC) $(FFLAGS) -o $@ $<

pymc_single_seg.o: pyMCsingleSegStime_NoLoop.f90
	$(FC) $(FFLAGS) -o $@ $<

reservoir:  precis.mod module_levelpool.o pymodule_levelpool.o

precis.mod: varPrecision.f90
	$(FC) $(FFLAGS) $<

install:
	cp *.o ../../../../src/python_routing_v02/fast_reach

%.o: %.F
	$(FC) $(FFLAGS) -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
