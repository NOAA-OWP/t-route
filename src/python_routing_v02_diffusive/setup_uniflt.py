from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system


# compile the fortran modules without linking
fortran_mod_comp =  'gfortran uniflowtzlt.f90 -c -o uniflowtzlt.o -O3 -fPIC'
system(fortran_mod_comp)

shared_obj_comp = 'gfortran pyuniflowtzlt.f90 -c -o pyuniflowtzlt.o -O3 -fPIC'
#print(f"shared_obj_comp: {shared_obj_comp}")
system(shared_obj_comp)

ext_modules = [Extension(# module name:
                         #"pygfunc",
                         "pyuniflowtzlt",                         
                         # source file:
                         #["pygfunc.pyx"],
                         ["pyuniflowtzlt.pyx"],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
			 #extra_compile_args=["-g"],
                         # other files to link to
                         extra_link_args=["uniflowtzlt.o", "pyuniflowtzlt.o"])]

setup(name = 'pyuniflowtzlt',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
