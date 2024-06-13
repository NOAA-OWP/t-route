# t-route - Tree-Based Channel Routing 

**Fast, flexible, modular channel routing for the National water model and beyond**:  

t-route, a dynamic channel routing model, offers a comprehensive solution for river network routing problems. It is designed to handle 1-D channel routing challenges in vector-based river network data, such as the USGS's NHDPlus High Resolution dataset, and OGC WaterML 2.0 Surface Hydrology Features (HY_Features) data model used in NextGen framework.

Provided a series lateral inflows for each node in a channel network, t-route computes the resulting streamflows. t-route requires that all routing computations srictly obey an upstream-to-downstream ordering. Such ordering facilitates the heterogenous application of routing models in a single river network. For example, hydrologic models - such as Muskingum Cunge - may be more appropriate for low-order headwater streams where backwater conditions have minimal impact on flooding. In contrast, t-route empowers users to apply computationally intensive hydraulic models, such as the diffusive wave approximation of St. Venant equations, for high-order streams and rivers, where backwater flooding is a significant concern. This flexibility enables users to allocate computational resources precisely where they are needed most, optimizing efficiency.

Expanding its capabilities, t-route now supports the OGC WaterML 2.0 Surface Hydrology Features (HY_Features) data model, facilitating the management of complex acyclic network connectivity. HY_Features data model provides users with a choice of routing solutions, including Muskingum-Cunge, Diffusive Wave, or Dynamic Wave.

**Key Features of t-route:**
- Adaptable to multiple network formulations.
- Supports a range of routing solutions.
- Redesigned core operations accessible via standard BMI functions.
- Modular design for independent reservoir and data assimilation modules.
- Capable of disaggregating routing across large spatial domains.
- Facilitates heterogenous application of routing models in a single network.

## Routing Model Comparisons
| Feature | Muskingum-Cunge (MC) | Diffusive Wave |
|---------|----------------------|----------------|
| Domain | Individual stream segment | Independent sub-network |
| Computing Method | Parallel on CONUS using junction order and short-time step | Serial on a sub-network using junction order |
| Routing Direction | Upstream to downstream | Both directions |
| Numerical Scheme | [One-dimensional explicit scheme](https://ral.ucar.edu/sites/default/files/public/WRF-HydroV5TechnicalDescription.pdf) | [Implicit Crank-Nicolson scheme](https://onlinelibrary.wiley.com/doi/10.1111/1752-1688.13080) |

t-route's flexible design is ideal for NOAA's National Water Model (NWM) 3.0, but its utility extends to various research and practical applications. While the code is currently bespoke for NWM, efforts are ongoing to generalize the core utilities, making t-route compatible with any standardized network and forcing data.

## General Sheme:
The figure below illustrates the workflow for executing t-route via two interfaces: the Command Line Interface (CLI) and the Basic Model Interface (BMI). When using the CLI, the user can select from two river network representations: NHDNetwork or HYFeatures. In contrast, the BMI exclusively supports HYFeatures. For the routing method, users have the option to apply either the Muskingum-Cunge method or the Diffusive Wave method.
<img src=https://raw.githubusercontent.com/NOAA-OWP/t-route/master/doc/images/scheme.png height=400>

## Project Overview:
- **Technology Stack**: Combines Python with Fortran for core routing model engines. The river network pre-processor, river network traversal framework, and time series data model are all written in python. The routing model engines (e.g. Muskingum-Cunge and diffusive wave) are primarily written in fortran, though we can imagine future additions to the engine suite being writted in any number of laguages that can be packaged as python extensions.
- **Current Status**: Focused on integration with NWM 3.0.
- **Demonstrations**: The `/test` directory includes a t-route demonstration on the Lower Colorado River, TX, showcasing capabilities like streamflow data assimilation and diffusive wave routing.

## Mission Alignment
t-route development is rigorously aligned with and guided by the NOAA Office of Water Prediction mission: *Collaboratively research, develop and deliver timely and consistent, state-of-the-science national hydrologic analyses, forecast information, data, guidance and equitable decision-support services to inform essential emergency management and water resources decisions across all time scales.*

## Structure of Folders:
- `src/kernel/`: Fortran modules.
- `src/troute-config/`: Configuration parser for t-route specific file.
- `src/troute-network/`: Manages network, data assimilation, and routing types.
- `src/troute-nwm/`: Coordinates t-route’s modules.
- `src/troute-routing/`: Handles flow segment and reservoir routing modules.

## Summary:
t-route represents streamflow channel routing and reservoir routing, assimilating data on vector-based channel networks. It fits into a broader framework where it interacts with land surface models, Forcing Engines, and coastal models, each playing a pivotal role in hydrologic forecasting and analysis.

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
geopandas
pyarrow
deprecated
```

## Installation

please see usage and testing below. Standby for docker container instructions in the near future.

## Configuration

Currently, there are no specific configuration details. Stand by for updates.

## Usage and Testing
To get a sense of the operation of the routing scheme, follow this sequence of commands:

```
# install required python modules
pip3 install numpy pandas xarray netcdf4 joblib toolz pyyaml Cython>3,!=3.0.4 geopandas pyarrow deprecated wheel

# clone t-toute
$ git clone --progress --single-branch --branch master http://github.com/NOAA-OWP/t-route.git

# compile and install
$ ./compiler.sh

# execute a demonstration test with NHD network
$ cd test/LowerColorado_TX
$ python3 -m nwm_routing -f -V4 test_AnA_V4_NHD.yaml

# OR

# execute a demonstration test with HYFeature network
$ cd test/LowerColorado_TX_v4
$ python3 -m nwm_routing -f -V4 test_AnA_V4_HYFeature.yaml
```

### T-Route Setup Instructions and Troubleshooting Guide for Windows Users

**Note**: The following instructions are for setting up T-Route on a Linux environment (standalone, no MPI). If you are using Windows, please install WSL (Windows Subsystem for Linux) before proceeding.

1. **Install Required Components:**
   - Open the WSL terminal.
   - Install Miniconda, Python, Pip, and Git.
   - Clone the Miniconda template repository: `git clone https://github.com/jameshalgren/miniconda-template.git`.
   - Follow the repository instructions to create a new environment.

2. **Activate the Environment:**
   - Activate the new environment created in the previous step.

3. **Install T-Route:**
   - Follow the instructions in the T-Route repository (https://github.com/NOAA-OWP/t-route/tree/master) for installation.
   - Ensure gcc, gfortran, and all required Python libraries are installed.

4. **NetCDF Issues:**
   - Resolve errors with NetCDF libraries (e.g., "netcdf.mod" not found) by running: `apt-get install *netcdf*`.
   - Locate the installed netcdf.mod (e.g., `find /usr/ -name *mod`).
   - Define the NETCDF path in the compiler.sh file in t-route (before the ‘if [-z “NETCDF …” ]’ statement): `export NETCDF="PATH_TO_NETCDF.MOD"`.

5. **Python Version:**
   - Define `alias python=python3` in the .bashrc file if `python` is not defined.
   - Change all instances of "python" to "python3" in the compiler file.

6. **Handle Permission Errors:**
   - Try compiling T-Route again after editing the compiler.
   - Use `sudo chmod 777 <path>` for "permission denied" errors (replace `<path>` with the relevant directory).

By following these instructions, you should successfully install and set up T-Route on your Linux system. For any issues or questions, feel free to seek assistance.


## Known issues

We are constantly looking to improve. Please see the Git Issues for additional information.

## Getting help

t-route team and the technical maintainers of the repository for more questions are: 
dongha.kim@noaa.gov 
sean.horvath@noaa.gov
amin.torabi@noaa.gov
jurgen.zach@noaa.gov

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
