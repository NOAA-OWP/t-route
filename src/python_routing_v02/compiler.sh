#TODO add options for clean/noclean, make/nomake, cython/nocython
#TODO include instuctions on blowing away entire package for fresh install e.g. rm -r ~/venvs/mesh/lib/python3.6/site-packages/troute/*
#set root folder of github repo (should be named t-route)
cd ../../
REPOROOT=`pwd`

#building reach and resevoir kernel files .o
cd $REPOROOT/src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS/
make clean
make
make install

#creates troute package
cd $REPOROOT/src/python_framework_v02
rm -rf build
#python setup.py --use-cython develop
python setup.py build_ext --inplace --use-cython
pip install -e .

#updates troute package with the execution script
cd $REPOROOT/src/python_routing_v02
rm -rf build
#python setup.py --use-cython develop
python setup.py build_ext --inplace --use-cython
pip install -e .
