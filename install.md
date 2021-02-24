# Installation instructions

This file provides more detailed instructions on how to install, configure, and get the project running. Summary instructions are provided in readme.md.

Step 1 -- Activate a python environment on your local machine 

Step 2 - Clone the repository 
``` 
git clone git@github.com:NOAA-OWP/t-route.git 
```

Step 3 - Compile necessary files
``` 
cd t-route/src/python_routing_v02/
./compiler.sh 
```

Step 4 - Run the model 

```
python3 compute_nhd_routing_SingleSeg_v02.py 
```

Step 4 - Alternative v01 model (defaults can be overwritten with custom input yaml)
``` 
cd ..
cd python_routing_v01
python3 compute_nhd_routing_SingleSeg_v01.py 
```
Examples:
```
python3 compute_nhd_routing_SingleSeg_v01.py -v -t
```
or
```
/t-route/test/input/yaml/CustomInput_short.yaml 
```
# Troubleshooting
- If you receive the following error.
    * FileNotFoundError: [Errno 2] No such file or directory: b'../../t-route/test/input/geo/Channels/RouteLink_NHDPLUS.nwm.v2.0.4.nc'

- Run the following command in ../../t-route/test/input/geo/Channels/
``` 
wget https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/nwm.v2.0.4/parm/domain/RouteLink_NHDPLUS.nc 
```

- Rename the downloaded file RouteLink_NHDPUS.nc to RouteLink_NHDPLUS.nwm.v2.0.4.nc

