#TODO add options for clean/noclean, make/nomake, cython/nocython
#TODO include instuctions on blowing away entire package for fresh install e.g. rm -r ~/venvs/mesh/lib/python3.6/site-packages/troute/*
#set root folder of github repo (should be named t-route)
cd ../../
REPOROOT=`pwd`
#For each build step, you can set these to true to make it build
#or set it to anything else (or unset) to skip that step
build_mc_kernel=true
build_diffusive_cnx_kernel=true
build_diffusive_cnt_kernel=true
build_reservoir_kernel=true
build_framework=true
build_routing=true

if [ -z "$F90" ]
then
    export F90="gfortran"
    echo "using F90='gfortran'"
fi

#if you have custom static library paths, uncomment below and export them
#export LIBRARY_PATH=<paths>:$LIBRARY_PATH
#if you have custom dynamic library paths, uncomment below and export them
#export LD_LIBRARY_PATHS=<paths>:$LD_LIBRARY_PATHS
export NETCDFINC=/usr/include/openmpi-x86_64/


if  [[ "$build_mc_kernel" == true ]]; then
  #building reach and resevoir kernel files .o
  cd $REPOROOT/src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS/
  make clean
  make  || exit
  make install || exit
fi

if  [[ "$build_diffusive_cnx_kernel" == true ]]; then
  #building reach and resevoir kernel files .o  
  cd $REPOROOT/src/fortran_routing/diffusive_v02frwk/
  make clean
  make diffusive.o
  make pydiffusive.o
  make install || exit
fi

if  [[ "$build_diffusive_cnt_kernel" == true ]]; then
  #building reach and resevoir kernel files .o  
  cd $REPOROOT/src/fortran_routing/diffusive_v02frwk/diffusive_cnt/
  make clean
  make arrays_module.o
  make matrix_module.o
  make var_module.o
  make arrays_section_module.o
  make xsec_attribute_module.o
  make constants_module.o
  make subtools.o
  make diffusive_cnt.o
  make pydiffusive_cnt.o
  make install || exit
fi

if [[ "$build_reservoir_kernel" == true ]]; then
  cd $REPOROOT/src/fortran_routing/Reservoir_Binding/Reservoirs
  make clean
  #make NETCDFINC=`nc-config --includedir` || exit
  #make binding_lp.a
  #make install_lp || exit
  make
  make install || exit

fi

if [[ "$build_framework" == true ]]; then
  #creates troute package
  cd $REPOROOT/src/python_framework_v02
  rm -rf build
  ##python setup.py --use-cython install
  ##python setup.py --use-cython develop
  python setup.py build_ext --inplace --use-cython || exit
  pip install -e . || exit
fi

if [[ "$build_routing" == true ]]; then
  #updates troute package with the execution script
  cd $REPOROOT/src/python_routing_v02
  rm -rf build
  #python setup.py --use-cython install
  #python setup.py --use-cython develop
  python setup.py build_ext --inplace --use-cython || exit
  pip install -e . || exit
fi
