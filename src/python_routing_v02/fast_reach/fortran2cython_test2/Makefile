# Modify these for your Fortran compiler:

# GFortran
F90 = gfortran
F90FLAGS = -Wall -std=f2008 -Wextra -Wimplicit-interface -fPIC
# Debug flags:
F90FLAGS += -g -fcheck=all
# Release flags:
#F90FLAGS += -O3 -march=native -ffast-math -funroll-loops

# Intel ifort
#F90 = ifort-12.0.191
#F90FLAGS = -stand f95 -warn all
# Debug flags:
#F90FLAGS += -check all
# Release flags:
#F90FLAGS += -xHOST -O3 -no-prec-div -static

#GCC
#/usr/include/python2.7
#/usr/lib64/python2.7/site-packages/numpy/core/include
#SPKG_LOCAL1=/usr
#SPKG_LOCAL2=/usr/lib64

#/usr/include/python3.6m
#/apd_common/channel_r/emanager/evaluation_manager_ir/eval_env_ir/lib/python3.6/site-packages/numpy/core/include
#SPKG_LOCAL1=/usr
#SPKG_LOCAL2=/apd_common/channel_r/emanager/evaluation_manager_ir/eval_env_ir

#/home/dongha.kim/anaconda3/include/python3.7m
#/home/dongha.kim/anaconda3/lib/python3.7/site-packages/numpy/core/include
SPKG_LOCAL1=/home/dongha.kim/anaconda3
SPKG_LOCAL2=/home/dongha.kim/anaconda3

CC = gcc
#CFLAGS = -Wall -fPIC -I$(SPKG_LOCAL1)/include/python2.7 -I$(SPKG_LOCAL2)/python2.7/site-packages/numpy/core/include
#CFLAGS = -Wall -fPIC -I$(SPKG_LOCAL1)/include/python3.6m -I$(SPKG_LOCAL2)/lib/python3.6/site-packages/numpy/core/include
CFLAGS = -Wall -fPIC -I$(SPKG_LOCAL1)/include/python3.7m -I$(SPKG_LOCAL2)/lib/python3.7/site-packages/numpy/core/include

fcython.so: fmesh_wrapper.o fcython.o fmesh.o utils.o types.o
	$(F90) -shared -o fcython.so $^ #-lpython3.7m

fcython.c: fcython.pyx
	cython fcython.pyx

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

fmesh_wrapper.o: fmesh.o
fmesh.o: types.o utils.o
utils.o: types.o
