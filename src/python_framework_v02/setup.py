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
    "troute.network.reach",
    sources=[
            "troute/network/reach.{}".format(ext),
            ],
    include_dirs=[np.get_include()],
    extra_objects=[],
    libraries=[],
    extra_compile_args=["-g"],
)

musk = Extension(
    "troute.network.musking.mc_reach",
    sources=[
            "troute/network/musking/mc_reach.{}".format(ext),
            ],
    include_dirs=[np.get_include(), "troute/network/"],
    extra_objects=[],
    libraries=[],
    extra_compile_args=["-g"],
)

reservoirs = Extension(
    "troute.network.reservoirs.levelpool.levelpool",
    sources=[
             "troute/network/reservoirs/levelpool/levelpool.{}".format(ext),
             ],
    include_dirs=[np.get_include(),  "troute/network/"],
    extra_objects=["./libs/binding_lp.a"],
    libraries=["gfortran", "netcdff", "netcdf"],
    extra_compile_args=["-g"],
)

package_data = {"troute.network": ["reach.pxd", "__init__.pxd"],
                "troute.network.musking": ["mc_reach.pxd", "__init__.pxd"],
                "troute.network.reservoirs":["__init__.pxd"],
                "troute.network.reservoirs.levelpool":["__init__.pxd", "levelpool.pxd"],
                 }
ext_modules = [reach, reservoirs, musk]

if USE_CYTHON:
    from Cython.Build import cythonize

    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": 3})

setup(
    name="troute.network",
    namespace_packages=["troute"],
    packages=find_namespace_packages(include=["troute.*"]),#["troute.network", "troute.network.reservoirs", "troute.network.reservoirs.levelpool"],
    ext_modules=ext_modules,
    ext_package="",
    package_data=package_data,
)
