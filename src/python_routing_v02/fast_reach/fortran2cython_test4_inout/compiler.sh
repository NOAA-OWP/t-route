cd ../../../../
REPOROOT=`pwd`
build_diff_kernel=true
build_routing=true

#if you have custom static library paths, uncomment below and export them
#export LIBRARY_PATH=<paths>:$LIBRARY_PATH
#if you have custom dynamic library paths, uncomment below and export them
#export LD_LIBRARY_PATHS=<paths>:$LD_LIBRARY_PATHS


if  [[ "$build_diff_kernel" == true ]]; then
  #building reach and resevoir kernel files .o  
  cd $REPOROOT/src/python_routing_v02/fast_reach/fortran2cython_test4_inout/fortran_code/
  make clean
  make diffusive.o
  make pydiffusive.o
  make install || exit
fi

if [[ "$build_routing" == true ]]; then
  #updates troute package with the execution script
  cd $REPOROOT/src/python_routing_v02/fast_reach/fortran2cython_test4_inout/
  rm -rf build
  #python setup.py --use-cython install
  #python setup.py --use-cython develop
  python setup.py build_ext --inplace --use-cython || exit
  pip install -e . || exit
fi