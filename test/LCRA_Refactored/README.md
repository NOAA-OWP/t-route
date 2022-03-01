# Lower Colorado AnA Test
This directory contains data from an operational 28-hour Natinal Water Model Analysis and Assimilation (AnA) run. The data have been subset to only include stream segments, waterbodies and gages within the Lower Colorado River basin, Texas, USA. This test demonstrates several key scientific capabilities of t-route, such as:
- Muskingum-Cunge channel routing
- Levelpool reseervoir routing
- Streamflow data assimilation
- RFC reservoir forecasting
- Muskingum-Cunge/Diffusive hybrid routing
The only capability not demonstrated by this test is USGS/USACE "hybrid" reservoir data assimilation. This is becuase the Lower Colorado basin does not contain any reservoirs for which USGS or USACE hybrid data assimiltion occurs. 

## Run the test
Simulation parameters for the test are contained in `t-route/test/LowerColorado_TX/test_AnA/yaml`. To execute the test:

```
bash

# navigate to python_routing_v02 directory
$ cd t-route/src/python_routing_v02

# compile routing kernels, reservoir modules, cython scripts, etc.
$ ./compiler.sh

# run the test
$ cd t-route/test/LowerColorado_TX
$ python -m nwm_routing -f -V3 test_AnA.yaml

```