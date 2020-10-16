# Test v02 with new edits to M-C routing module
## What changed?
- At the core of this pull request are changes to `MCsingleSegStime_f2py_NOLOOP.f90`. Now this `muskingcungenwm` subroutine in this module returns three diagnostic variables: kinematic celerity, courant number, and M-C X parameter. 
- To accommodate this change with v02 routing, I edited `pyMCsingleSegStime_NoLoop.f90`, to accept three additional output variables from `muskingcungenwm`, while only returning qdc, velc, and depthc. As such, the `mc_reach.compute_network` and specifically `mc_reach.compute_reach_kernel` will run just fine. However, while additional diagnostic variables provided by `muskingcungenwm`, they are not yet "seen" or made available to the `mc_reach`.
- To allow testing of v02, I also edited `compute_nhd_routing_SingleSeg_v02.py` to add the `/fast_reach` directory to the system path, else it is unable to find and load `mc_reach`.

## How to test that Adam did not break v02
### 0. Clone this branch, or figure out how to work in Google Colab.
### 1. Compile v02 cython
- While the compile commands will be unique to your system requirements, here is what works with A. Wlostowski's system (a shell script):

```
# define directory path variables, assumed root = /t-route
FORTRAN_ROUTING_DIR="/src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS"
FORTRAN_RES_DIR="/src/fortran_routing/mc_pylink_v00/Reservoir_singleTS"
FAST_REACH_DIR="/src/python_routing_v02/fast_reach"

cd $FORTRAN_ROUTING_DIR
gfortran varPrecision.f90 -c -O0 -fPIC
gfortran -c -O0 -fPIC -o mc_single_seg.o MCsingleSegStime_f2py_NOLOOP.f90
gfortran pyMCsingleSegStime_NoLoop.f90 -c -o pymc_single_seg.o -O3 -fPIC
cp $FORTRAN_ROUTING_DIR/*.o $FAST_REACH_DIR 

cd $FORTRAN_RES_DIR
gfortran varPrecision.f90 -c -O0 -fPIC
gfortran -c -O2 -fPIC -o module_levelpool.o module_levelpool.f90
gfortran -c -O2 -fPIC -o pymodule_levelpool.o pymodule_levelpool.f90
cp $FORTRAN_RES_DIR/*.o $FAST_REACH_DIR 

# define more directory path variables - system specific!
numpy_I="/home/user/tr/lib/python3.8/site-packages/numpy/core/include"
py_I="/usr/include/python3.8"
py_lib="/home/user/tr/lib"

cd $FAST_REACH_DIR
cython -3 -v -p --gdb --line-directives -Wextra --cleanup 3 fortran_wrappers.pxd *.pyx
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o fortran_wrappers.o fortran_wrappers.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o mc_reach.o mc_reach.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o reservoir.o reservoir.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o reach.o reach.c
gcc -pthread -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -fno-strict-aliasing -I$numpy_I -I$py_I -mtune=generic -c -o utils.o utils.c
gcc -pthread -shared -L $py_lib -lgfortran -o reach.cpython-38-x86_64-linux-gnu.so mc_single_seg.o pymc_single_seg.o fortran_wrappers.o reach.o
gcc -pthread -shared -L $py_lib -lgfortran -o reservoir.cpython-38-x86_64-linux-gnu.so module_levelpool.o pymodule_levelpool.o fortran_wrappers.o reservoir.o
gcc -pthread -shared -L $py_lib -o mc_reach.cpython-38-x86_64-linux-gnu.so mc_reach.o
gcc -pthread -shared -L $py_lib -o utils.cpython-38-x86_64-linux-gnu.so utils.o

```

### 2. Run the test 

```
cd src/python_routing_v02
python compute_nhd_routing_SingleSeg_v02.py
```