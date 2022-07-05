# t-route - Tree-Based Channel Routing 

**Fast, flexible, modular channel routing for the National water model and beyond**:  

t-route solves 1-D channel routing problems for vector-based river network data, such as USGS's NHDPlus High Resolution dataset. Provided a series lateral inflows for each node in a channel network, t-route computes the resulting streamflows. t-route requires that all routing computations srictly obey an upstream-to-downstream ordering. Such ordering facilitates the heterogenous application of routing models in a single river network. For example, hydrologic models - such as Muskingum Cunge - may be more appropriate for low-order headwater streams where backwater conditions are less likely to affect regional flooding. On the other hand, more computationally intensive hydraulic models - such as the diffusive wave approximation of St. Venant equations - are better suited for high-order streams and rivers where backwater flooding regularly impacts regional communities. With t-route, users are free to create routing models that use hydraulic solutions (and their computational resources) only where they are most needed, while relying on simpler hydrologic solutions elsewhere. 

t-route is designed for use in NOAA's National Water Model (NWM) 3.0, though can also be leveraged for research and practial applications outside of the national water model. At the moment, much of our code is bespoke to NWM application, and we are constantly working to generisize and abstract the core utilities of the codebase so that t-route can be used with any standarized network and forcing data.

t-route development is rigorously aligned with and guided by the NOAA Office of Water Prediction mission: *Collaboratively research, develop and deliver timely and consistent, state-of-the-science national hydrologic analyses, forecast information, data, guidance and equitable decision-support services to inform essential emergency management and water resources decisions across all time scales.*  

The project and program contain the following elements. 
  - **Technology stack**: The river network pre-processor, river network traversal framework, and time series data model are all written in python. The routing model engines (e.g. Muskingum-Cunge and diffusive wave) are primarily written in fortran, though we can imagine future additions to the engine suite being writted in any number of laguages that can be packaged as python extensions. 
  - **Status**:  The project is currently focussed on preparing t-route for application NWM 3.0.
  - **Demos**: The `/test` directory contains a canned t-route demonstration of routing on the Lower Colorado River, TX for a standard NWM Analsis and Assimilation simulation. A suite of operationally useful capabilities are demonstrated, including streamflow data assimilation, diffusive wave routing along the Lower Colorado mainstem, waterbody simulation, and RFC reservoir forecasting. 

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
