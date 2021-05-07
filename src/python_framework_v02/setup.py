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
    include_dirs=[np.get_include()],
    extra_objects=[],
    libraries=[],
    extra_compile_args=["-g"],
)

levelpool_reservoirs = Extension(
    "troute.network.reservoirs.levelpool.levelpool",
    sources=[
             "troute/network/reservoirs/levelpool/levelpool.{}".format(ext),
             ],
    include_dirs=[np.get_include()],
    extra_objects=["./libs/binding_lp.a"],
    libraries=["gfortran", "netcdff", "netcdf"],
    extra_compile_args=["-g"],
)

hybrid_reservoirs = Extension(
    "troute.network.reservoirs.hybrid.hybrid",
    sources=[
             "troute/network/reservoirs/hybrid/hybrid.{}".format(ext),
             ],
    include_dirs=[np.get_include()],
    extra_objects=["./libs/bind_hybrid.a"],
    libraries=["gfortran", "netcdff", "netcdf"],
    extra_compile_args=["-g"],
)

rfc_reservoirs = Extension(
    "troute.network.reservoirs.rfc.rfc",
    sources=[
             "troute/network/reservoirs/rfc/rfc.{}".format(ext),
             ],
    include_dirs=[np.get_include()],
    extra_objects=["./libs/bind_rfc.a"],
    libraries=["gfortran", "netcdff", "netcdf"],
    extra_compile_args=["-g"],
)

package_data = {"troute.network": ["reach.pxd", "__init__.pxd"],
                "troute.network.musking": ["mc_reach.pxd", "__init__.pxd"],
                "troute.network.reservoirs":["__init__.pxd"],
                "troute.network.reservoirs.levelpool":["__init__.pxd", "levelpool.pxd"],
                "troute.network.reservoirs.hybrid":["__init__.pxd", "hybrid.pxd"],
                "troute.network.reservoirs.rfc":["__init__.pxd", "rfc.pxd"],
                 }
ext_modules = [reach, levelpool_reservoirs, hybrid_reservoirs, rfc_reservoirs, musk]

if USE_CYTHON:
    from Cython.Build import cythonize

    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": 3})

setup(
    name="Network",
    packages=["troute", "troute.network", "troute.network.reservoirs", "troute.network.reservoirs.levelpool", "troute.network.reservoirs.hybrid", "troute.network.reservoirs.rfc"],
    ext_modules=ext_modules,
    ext_package="",
    package_data=package_data,
)
