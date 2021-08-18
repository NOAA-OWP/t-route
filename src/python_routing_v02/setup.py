from setuptools import setup, find_namespace_packages
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

reach = Extension(
    "troute.routing.fast_reach.reach",
    sources=["troute/routing/fast_reach/reach.{}".format(ext)],
    extra_objects=[
        "troute/routing/fast_reach/mc_single_seg.o",
        "troute/routing/fast_reach/pymc_single_seg.o",
    ],
    extra_compile_args=["-g"],
)

mc_reach = Extension(
    "troute.routing.fast_reach.mc_reach",
    sources=["troute/routing/fast_reach/mc_reach.{}".format(ext)],
    include_dirs=[np.get_include()],
    libraries=[],
    library_dirs=[],
    extra_objects=[],
    extra_compile_args=["-g"],
)

diffusive = Extension(
    "troute.routing.fast_reach.diffusive",
    sources=["troute/routing/fast_reach/diffusive.{}".format(ext)],
    extra_objects=[
        "troute/routing/fast_reach/diffusive.o",
        "troute/routing/fast_reach/pydiffusive.o",
    ],
    include_dirs=[np.get_include()],
    extra_compile_args=["-g"],
    libraries=["gfortran"],
)

diffusive_cnt = Extension(
    "troute.routing.fast_reach.diffusive_cnt",
    sources=["troute/routing/fast_reach/diffusive_cnt.{}".format(ext)],
    extra_objects=[
        "troute/routing/fast_reach/arrays_module.o",
        "troute/routing/fast_reach/matrix_module.o",
        "troute/routing/fast_reach/var_module.o",
        "troute/routing/fast_reach/arrays_section_module.o",
        "troute/routing/fast_reach/xsec_attribute_module.o",
        "troute/routing/fast_reach/constants_module.o",
        "troute/routing/fast_reach/subtools.o",
        "troute/routing/fast_reach/diffusive_cnt.o",
        "troute/routing/fast_reach/pydiffusive_cnt.o",
    ],
    include_dirs=[np.get_include()],
    extra_compile_args=["-g"],
    libraries=["gfortran"],
)

package_data = {"troute.fast_reach": ["reach.pxd", "fortran_wrappers.pxd", "utils.pxd"]}
ext_modules = [reach, mc_reach, diffusive, diffusive_cnt]

if USE_CYTHON:
    from Cython.Build import cythonize

    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": 3})

setup(
    name="troute.routing",
    namespace_packages=["troute"],
    packages=find_namespace_packages(
        include=["troute.*"]
    ),  # ["troute.fast_reach", "troute.routing"],
    ext_modules=ext_modules,
    ext_packages="",
    package_data=package_data,
)
