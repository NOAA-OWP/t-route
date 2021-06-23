# NWM Tree-based Inland Hydraulic Routing Project 

**Fast, flexible, modular channel routing for the National water model**:  


Scaling the hydraulic routing problem to the size of the U.S. National Water Model required a new approach. The existing routing algorithm does not allow for incorporation of downstream conditions in the calculation of stream flow and depth. A hydraulic model which explicitly considers these influences is required for simulation of regional flooding where waterways begin to influence one another as well as when a major flood causes backwater flooding into a tributary. The program under development here seeks to effectively manage the traversal of a network of streams with defined hydraulic properties specifically for the purpose of hydraulic routing in an operational flood and water resources forecasting system. The principles of graph development and traversal as they are applied (imperfectly) here are possibly applicable to problems in other areas. 

In particular, given the observation that the routing of independent river networks may be completely decoupled, we have worked to create a framework for the routing computation which uses a knowledge of the topological relationship of segments, reaches, and junctions within a river network to separate the computation into parallelizable portions, as shown in the following image. 

![](./doc/bluecyan.gif)

The project and program contain the following elements. 
  - **Technology stack**: The hydrofabric pre-processor, river network traversal framework, and time series data model are all written in python. The routing model engines are primarily written in fortran, but we are also experimenting with ML-based methods and can imagine any number of options called from within the network traversal framework.
<<<<<<< HEAD
  - **Status**:  The project is currently in development phase 2 -- we are making the first connections of the various components within a single module. Phase 3 shoulud begin around July 2020 and we will be working demonstrations of the framework with operational outputs bootstrapped from the current national water model. Eventually, there will be a [CHANGELOG](CHANGELOG.md).
=======
  - **Status**:  The project is currently in development phase 2 -- we are making the first connections of the various components within a single module. Phase 3 shoulud begin around July 2020 and we will be working demonstrations of the framework with operational outputs bootstrapped from the current national water model. For additional detail on versions, please see the [CHANGELOG](changelog.md).
>>>>>>> upstream/master
  - **Demos**: The `notebooks` folder has a number of python notebooks, many of which can be executed from within the Google colaboratory environment, which demonstrate various aspects of the project. 

In general, as with all things related to the National Water Center, our goal is to save
more lives and prevent more property damage by more accurately forecasting the state of the nations
rivers and other water resources, especially with respect to floods and by getting that
information more quickly and understandably to people who can and will do something 
useful with it to those ends.

The purpose in sharing this via GitHub is to efficiently make possible collaborative 
editing -- because there are certainly those who have more experience and expertise who
might comment on how this has been done or how it might be better done here. 

## Configuration and Dependencies

This program uses the following system packages:
```
python3
python3-devel
gcc-gfortran
```

<<<<<<< HEAD
... and the following python modules:
=======
... and the following non-default python modules:
>>>>>>> upstream/master
```
geopandas 
numpy 
pandas 
xarray 
netcdf4 
<<<<<<< HEAD
=======
tqdm # for a download progress bar when pulling test datasets from the NOAA site.
>>>>>>> upstream/master
```

## Installation

please see usage and testing below. Standby for docker container instructions in the near future.

## Configuration

Presently, there are no specific configuration details. Stand by.

## Usage and Testing
The following sequence of commands should provide a sense of the operation of the routing scheme:

```
pip3 install geopandas numpy pandas xarray netcdf4 
git clone --progress --single-branch --branch master http://github.com/NOAA-OWP/t-route.git
cd src/fortran_routing/mc_pylink_v00/MC_singleCH_singleTS/
f2py3 -c varSingleChStime_f2py.f90  MCsingleChStime_f2py_clean.f90  -m mc_sc_stime
cd -
cd src/python_routing_v02/
# Execute a serial (~2.5 minutes) and a parallel test (~0.5 minutes)
# times from 6 cores (x2 threads per core), 3.7GHz
<<<<<<< HEAD
python3 compute_nhd_routing.py; python3 parallel_compute_nhd_routing.py.
=======
python3 compute_nhd_routing.py; python3 parallel_compute_nhd_routing.py
>>>>>>> upstream/master
```
 
In addition, the notebooks in the `notebooks` folder represent experiments to demonstrate the base classes and basic functionality of libraries contained in rest of the repoository. Collectively they are intended to support incorporation of evaluation of the dynamic routing codes and related utilities being considered and prepared as potential candidates for use as the base dynamic hydraulic engine for the NWM. 

Most of the notebooks may be opened in either Jupyter or with Google Colaboratory.

## Known issues

We are constantly looking to improve. Please see the Git Issues for additional information.

## Getting help

If you have any questions, please contact james.halgren@noaa.gov or dongha.kim@noaa.gov, the technical maintainers of the repository. 

## Getting involved

<<<<<<< HEAD
Our current focus is improving the speed of the parallel tree traversal. We welcome your thoughts, recommendations, comments, and of course, PRs. 

Please feel free to fork the repository and let us know if you will be issuing a pull request. 
More instructions will eventually be documented in [CONTRIBUTING](CONTRIBUTING.md).
=======
Among other things, we are working on preparing more robust I/O capability for testing and on improving the speed of the parallel tree traversal. We welcome your thoughts, recommendations, comments, and of course, PRs. 

Please feel free to fork the repository and let us know if you will be issuing a pull request. 
More instructions will eventually be documented in [CONTRIBUTING](contributing.md).
>>>>>>> upstream/master


----

## Open source licensing info
1. [TERMS](TERMS.md)
2. [LICENSE](LICENSE)


----

## Credits and references

A great deal of credit is owed to Drs. Ehab Meslehe and Fred Ogden, and as well to the entire NCAR WRF-Hydro development team.

----
