#How to build and execute the cython-based routing function

The following are the dependencies required for the cython-based call. 
`pip install cython numpy joblib`

The following are needed for the basic data io. 
`pip install geopandas xarray tqdm requests netcdf4`
(geopandas installs pandas)

The python development headers are needed if they are not installed. The scripts below assume that the headers are linked to $VIRTUAL_ENV/include/<python_version>
`sudo yum install python3-devel`
(or appropriate substitute package manager)

The `$VIRTUAL_ENV` variables below can be replaced by `$CONDA_ENV` in most cases and will work fine. 

The python version and other system variables are encoded in the include folders and in the shared object name below. This short python script accesses that information. 
```
import sysconfig
sysconfig.get_config_var('SOABI')
# prints e.g., 'cpython-36m-x86_64-linux-gnu'
```

Compilation with gfortran and gcc
```
cd src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS/
gfortran varPrecision.f90 -c -o var_precision.o -O3 -fPIC
gfortran MCsingleSegStime_f2py_NOLOOP.f90 -c -o mc_single_seg.o -O3 -fPIC
gfortran pyMCsingleSegStime_NoLoop.f90 -c -o pymc_single_seg.o -O3 -fPIC
cp *.o ../../../../src/python_routing_v02
cd ../../../../src/python_routing_v02
cython -3 -v -p --line-directives -Wextra --cleanup 3 mc_reach.pyx
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -g -fwrapv -O3 -fno-strict-aliasing -Wall -Wstrict-prototypes -fPIC -I$VIRTUAL_ENV/lib/python3.6/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.6m -c mc_reach.c -o mc_reach.o
gcc -pthread -shared -L$VIRTUAL_ENV/lib var_precision.o mc_single_seg.o pymc_single_seg.o mc_reach.o -lgfortran -o mc_reach.cpython-36m-x86_64-linux-gnu.so
python3 compute_nhd_routing_SingleSeg_arr2.py
```


For compilation with ifort and icc (such as on the Cheyenne supercomputer)
```
cd src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS/
ifort varPrecision.f90 -c -o var_precision.o -O3 -fPIC
ifort MCsingleSegStime_f2py_NOLOOP.f90 -c -o mc_single_seg.o -O3 -fPIC
ifort pyMCsingleSegStime_NoLoop.f90 -c -o pymc_single_seg.o -O3 -fPIC
cp *.o ../../../../src/python_routing_v02
cd ../../../../src/python_routing_v02
cython -3 -v -p --line-directives -Wextra --cleanup 3 mc_reach.pyx
icc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I$VIRTUAL_ENV/lib/python3.7/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.7m -c mc_reach.c -o mc_reach.o
icc -pthread -shared -L$VIRTUAL_ENV/lib var_precision.o mc_single_seg.o pymc_single_seg.o mc_reach.o -Iifort -o mc_reach.cpython-37m-x86_64-linux-gnu.so
python3 compute_nhd_routing_SingleSeg_arr2.py
```
