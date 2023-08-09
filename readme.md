# t-route - Tree-Based Channel Routing 

**Fast, flexible, modular channel routing for the National water model and beyond**:  

t-route solves 1-D channel routing problems for vector-based river network data, such as USGS's NHDPlus High Resolution dataset. Provided a series lateral inflows for each node in a channel network, t-route computes the resulting streamflows. t-route requires that all routing computations srictly obey an upstream-to-downstream ordering. Such ordering facilitates the heterogenous application of routing models in a single river network. For example, hydrologic models - such as Muskingum Cunge - may be more appropriate for low-order headwater streams where backwater conditions are less likely to affect regional flooding. On the other hand, more computationally intensive hydraulic models - such as the diffusive wave approximation of St. Venant equations - are better suited for high-order streams and rivers where backwater flooding regularly impacts regional communities. With t-route, users are free to create routing models that use hydraulic solutions (and their computational resources) only where they are most needed, while relying on simpler hydrologic solutions elsewhere. 

t-route is designed for use in NOAA's National Water Model (NWM) 3.0, though can also be leveraged for research and practial applications outside of the national water model. At the moment, much of our code is bespoke to NWM application, and we are constantly working to generisize and abstract the core utilities of the codebase so that t-route can be used with any standarized network and forcing data.

t-route development is rigorously aligned with and guided by the NOAA Office of Water Prediction mission: *Collaboratively research, develop and deliver timely and consistent, state-of-the-science national hydrologic analyses, forecast information, data, guidance and equitable decision-support services to inform essential emergency management and water resources decisions across all time scales.*  

The project and program contain the following elements. 
  - **Technology stack**: The river network pre-processor, river network traversal framework, and time series data model are all written in python. The routing model engines (e.g. Muskingum-Cunge and diffusive wave) are primarily written in fortran, though we can imagine future additions to the engine suite being writted in any number of laguages that can be packaged as python extensions. 
  - **Status**:  The project is currently focussed on preparing t-route for application NWM 3.0.
  - **Demos**: The `/test` directory contains a canned t-route demonstration of routing on the Lower Colorado River, TX for a standard NWM Analsis and Assimilation simulation. A suite of operationally useful capabilities are demonstrated, including streamflow data assimilation, diffusive wave routing along the Lower Colorado mainstem, waterbody simulation, and RFC reservoir forecasting. 

##
##
## INSTRUCTIONS TO USE ON WSL (Windows Subsystem for Linux), SPECIFICALLY ON UBUNTU 20.04 (Long Term Stable distribution)
## ADDED 09-Aug-2023 (legacy instructions see below)
##
##

## Install Ubuntu on WSL2 (any Windows 10/11 PC; tested on Windows 11)

Open Windows PowerShell, command:

>> wsl --install Ubuntu-20.04

The Linux distribution will be installed and tracked by a progress bar in the Power Shell window. Once done, you will be prompted to enter a Linux username and password - please do so. The distribution is automatically, and you will be logged into a Linux terminal. 


For verification, you can open a new PowerShell, and verify creation of the distribution using the following command:

>> wsl -l -v

, which should result in something like the following (where the line with Ubuntu-20.04 should have been added; other lines are optional):

  NAME             STATE           VERSION
* rhel79           Stopped         2
  Ubuntu-20.04     Running         2
  Ubuntu_v2        Stopped         2
  Ubuntu           Stopped         2
  Debian           Running         2
  Ubuntu_Server    Stopped         2


##
## In Ubuntu: perform updates
##

>> sudo apt update

[OPTIONAL: If there are strange error messages about updates not being valid yet, that's due to the system clock, and can be fixed by synchronizing the hardware clock:
>> sudo hwclock --hctosys ]

>> sudo apt -y upgrade


##
## Install Anaconda (in Ubuntu)
##

Check the archive in  https://repo.continuum.io/archive for available releases. Comment: a stable older version is given as an example here.

Download the distribution:

>> wget https://repo.continuum.io/archive/Anaconda3-2022.05-Linux-x86_64.sh


Install Anaconda - Review and confirm license agreement. As installation directory, installed in home directory /home/jz/anaconda3 instead of accepting default installation directory (/root/anaconda3). When done,  agree to running conda init.

>> sudo bash ./Anaconda3-2023.07-1-Linux-x86_64.sh


Verify that the python installation is the Anaconda one, not the system one. Issuing the following command:
>> which python
should give something like /home/[YOUR LINUX USER NAME]/anaconda3/bin/python, NOT /use/bin/python or similar

The following commands' purpose is to be able to use the anaconda-navigator GUI:
>> sudo apt-get install libxss1
>> sudo apt install libegl1-mesa libegl1
>> sudo apt install x11-apps


##
## Create conda environment: use Python 3.9 (in line with existing installations)
##

Create a local conda  environment:
>> conda create -n py39 python=3.9

Activate it:
>> conda activate py39

Install t-route python packages using conda to make sure there are no conflicts (pip would be much faster, but the installation may not run):
>> conda install -c conda-forge numpy
>> conda install -c conda-forge pandas
>> conda install -c conda-forge xarray
>> conda install -c conda-forge netcdf4
>> conda install -c conda-forge joblib
>> conda install -c conda-forge toolz
>> conda install -c conda-forge Cython
>> conda install -c conda-forge pyyaml
>> conda install -c conda-forge geopandas
>> conda install -c conda-forge pyarrow
>> conda install -c conda-forge deprecated

Install local gcc and gfortran in conda environment:
>> conda install -c conda-forge gcc
>> conda install -c conda-forge gfortran

Prepare compilation:
>> conda install -c conda-forge make
>> conda install -c conda-forge netcdf-fortran


##
## It is good practice to deactivate and activate the environment again and run any .bashrc additions in the meantime:
##

>> conda deactivate
>> cd
>> source .bashrc
>> conda activate py39


##
## Clone the t-route Github repository
##

Recommended: create unique directory, cd into it:
>> mkdir tRoute_39_2
>> cd tRoute_39_2

Clone t-route:
>> git clone --progress --single-branch --branch master http://github.com/NOAA-OWP/t-route.git
>> cd t-route

Set some system variables (recommended to do that in, e.g., .bashrc if used more than a couple of times:
>> export NETCDFINC=/home/jz/.conda/envs/py39/include/
>> export WSL=1

The code should then compile when running:
>> ./compiler.sh




##
##
## LEGACY INSTRUCTIONS FOLLOW
##
##

## Configuration and Dependencies

This program uses the following system packages:
```
python3
gcc-gfortran
```

... and the following non-default python modules:
``` 
numpy 
pandas 
xarray 
netcdf4 
joblib
toolz
Cython
pyyaml
```

## Installation

please see usage and testing below. Standby for docker container instructions in the near future.

## Configuration

Presently, there are no specific configuration details. Stand by.

## Usage and Testing
The following sequence of commands should provide a sense of the operation of the routing scheme:

```
# install required python modules
$ pip3 install numpy pandas xarray netcdf4 joblib toolz pyyaml Cython

# clone t-toute
$ git clone --progress --single-branch --branch master http://github.com/NOAA-OWP/t-route.git

# compile and install
$ ./compiler.sh

# execute a demonstration test
$ cd test/LowerColorado_TX
$ python3 -m nwm_routing -f test_AnA.yaml
```

## Known issues

We are constantly looking to improve. Please see the Git Issues for additional information.

## Getting help

If you have any questions, please contact adam.wlostowski@noaa.gov or dongha.kim@noaa.gov, the technical maintainers of the repository. 

## Getting involved

Among other things, we are working on preparing more robust I/O capability for testing and on improving the speed of the parallel tree traversal. We welcome your thoughts, recommendations, comments, and of course, PRs. 

Please feel free to fork the repository and let us know if you will be issuing a pull request. 
More instructions will eventually be documented in [CONTRIBUTING](contributing.md).


----

## Open source licensing info
1. [TERMS](TERMS.md)
2. [LICENSE](LICENSE)


----

## Credits and references

A great deal of credit is owed to Drs. Ehab Meslehe and Fred Ogden, and as well to the entire NCAR WRF-Hydro development team. Continued leadership and support from Dr. Trey Flowers.

----
