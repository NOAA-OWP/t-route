from setuptools import setup
from distutils.extension import Extension
import sys
import numpy as np

"""
If source ships with the cython generated .c code, then cython isn't a hard requirement
This setup can use if if passed the --use-cython flag, otherwise it looks for the c source
and builds using distutils.
"""

if "--use-cython" in sys.argv:
    USE_CYTHON = True
    sys.argv.remove("--use-cython")
else:
    USE_CYTHON = False

ext = "pyx" if USE_CYTHON else "c"

diffusive = Extension(
    "diffusive",
    sources=["diffusive.{}".format(ext)],
    extra_objects=["diffusive.o", "pydiffusive.o"],
    extra_compile_args=["-g"],
    libraries = ['gfortran'],
)

ext_modules = [diffusive]

if USE_CYTHON:
    from Cython.Build import cythonize

    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": 3})

setup(
    name="compute_diffusive",
    ext_modules=ext_modules,
)
