from setuptools import setup
from distutils.extension import Extension
import sys
import numpy as np
"""
If source ships with the cython generated .c code, then cython isn't a hard requirement
This setup can use if if passed the --use-cython flag, otherwise it looks for the c source
and builds using distutils.
"""

if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False

ext = 'pyx' if USE_CYTHON else 'c'

reach = Extension("troute.network.reach",
        sources = ["troute/network/reach.{}".format(ext)],
        include_dirs=[np.get_include()],
        extra_objects = [],
        extra_compile_args = ['-g'],
        )

package_data = {
    'troute.network':['reach.pxd']
    }
ext_modules=[
    reach
    ]

if USE_CYTHON:
    from Cython.Build import cythonize
    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": 3})

setup(
  name = 'Network',
  packages=['troute', 'troute.network'],
  ext_modules = ext_modules,
  package_data = package_data,
)
