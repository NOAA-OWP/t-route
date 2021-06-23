#How to build and execute the cython-based routing function

The following are the dependencies required for the cython-based call. 
`pip install cython numpy joblib`

The following are needed for the basic data io. 
`pip install geopandas xarray tqdm requests netcdf4`
(geopandas installs pandas)

The python development headers are needed if they are not installed. The scripts below assume that the headers are linked to $VIRTUAL_ENV/include/<python_version>
`sudo yum install python3-devel`
(or appropriate substitute package manager)

The `$VIRTUAL_ENV` variables below can be replaced by `$CONDA_PREFIX` in most cases and will work fine. 

The python version and other system variables are encoded in the include folders and in the shared object name below. This short python script accesses that information. 
```
import sysconfig
sysconfig.get_config_var('SOABI')
# prints e.g., 'cpython-36m-x86_64-linux-gnu'
```

A shell script with the following commands is sufficient for compilation with gfortran and gcc, and the command `python3 compute_nhd_routing_SingleSeg_arr2.py` will demonstrate the application. 
```
set -x

gfortran varPrecision.f90 -c -O0 -fPIC
gfortran -c -O2 -fPIC -o mc_single_seg.o MCsingleSegStime_f2py_NOLOOP.f90
gfortran -c -O2 -fPIC -o pymc_single_seg.o pyMCsingleSegStime_NoLoop.f90
cp *.o ../../../../src/python_routing_v02/fast_reach
cd ../Reservoir_singleTS/
gfortran -c -O2 -fPIC -o module_levelpool.o module_levelpool.f90
gfortran -c -O2 -fPIC -o pymodule_levelpool.o pymodule_levelpool.f90
cp *.o ../../../../src/python_routing_v02/fast_reach
cd ../../../../src/python_routing_v02/fast_reach


cython -3 -v -p --gdb --line-directives -Wextra --cleanup 3 fortran_wrappers.pxd *.pyx

gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$VIRTUAL_ENV/lib/python3.6/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.6m -mtune=generic -c -o fortran_wrappers.o fortran_wrappers.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$VIRTUAL_ENV/lib/python3.6/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.6m -mtune=generic -c -o mc_reach.o mc_reach.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$VIRTUAL_ENV/lib/python3.6/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.6m -mtune=generic -c -o reservoir.o reservoir.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$VIRTUAL_ENV/lib/python3.6/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.6m -mtune=generic -c -o reach.o reach.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$VIRTUAL_ENV/lib/python3.6/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.6m -mtune=generic -c -o utils.o utils.c

gcc -pthread -shared -L$VIRTUAL_ENV/lib -lgfortran -o reach.cpython-36m-x86_64-linux-gnu.so mc_single_seg.o pymc_single_seg.o fortran_wrappers.o reach.o
gcc -pthread -shared -L$VIRTUAL_ENV/lib -lgfortran -o reservoir.cpython-36m-x86_64-linux-gnu.so module_levelpool.o pymodule_levelpool.o fortran_wrappers.o reservoir.o
gcc -pthread -shared -L$VIRTUAL_ENV/lib -o mc_reach.cpython-36m-x86_64-linux-gnu.so mc_reach.o
gcc -pthread -shared -L$VIRTUAL_ENV/lib -o utils.cpython-36m-x86_64-linux-gnu.so utils.o
cp *.so ../
```


For compilation with ifort and icc (such as on the Cheyenne supercomputer), `ifort` may be substituted for `gfortran`. 
The compiling and linking steps for example are as shown:
```
icc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$VIRTUAL_ENV/lib/python3.7/site-packages/numpy/core/include -I$VIRTUAL_ENV/include/python3.7m -c mc_reach.c -o mc_reach.o
icc -pthread -shared -L$VIRTUAL_ENV/lib mc_reach.o -Iifort -o mc_reach.cpython-37m-x86_64-linux-gnu.so
```
