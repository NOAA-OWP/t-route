from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


#Document to set this up
#When installing this module from source
#$LD_LIBRARY_PATH=/path/to/my/netcdf/:$LD_LIBRARY_PATH setup.py


ext_modules=[
    Extension("compute_kernel_lp", 
              sources = ["compute_kernel_lp.pyx"],
              include_dirs = [],
              libraries=['gfortran', 'netcdff', 'netcdf'],
              library_dirs=['/lib/gcc/x86_64-redhat-linux/4.8.5', '/usr/lib64/openmpi/lib/'],
              extra_objects=['fortran/Reservoirs/binding.a']),

    Extension("compute_kernel_hybrid", 
              sources = ["compute_kernel_hybrid.pyx"],
              include_dirs = [],
              libraries=['gfortran', 'netcdff', 'netcdf'],
              library_dirs=['/lib/gcc/x86_64-redhat-linux/4.8.5', '/usr/lib64/openmpi/lib/'],
              extra_objects=['fortran/Reservoirs/bind_hybrid.a']),

    Extension("compute_kernel_rfc", 
              sources = ["compute_kernel_rfc.pyx"],
              include_dirs = [],
              libraries=['gfortran', 'netcdff', 'netcdf'],
              library_dirs=['/lib/gcc/x86_64-redhat-linux/4.8.5', '/usr/lib64/openmpi/lib/'],
              extra_objects=['fortran/Reservoirs/bind_rfc.a'])
]

setup(
  name = 'RoutingKernels',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
)
