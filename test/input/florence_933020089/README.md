# Hurricane Florence Test Case
This directory contains data from a WRF-Hydro simulation
of Hurricane Florence on a domain in North Carolina, USA.
Because WRF-Hydro output data can be quite bulky, the
lateral inflow and streamflow data found here have been
constructed by pre-processing, then re-exporting the
contents of the CHRTOUT channel output and HYDRO_RST
WRF-Hydro restart files as follows:
- *.CHRTOUT files. -->
    test/input/florence_933020089/OUTPUT/SeptemberFlorenceQlateral_wDates.csv
    test/input/florence_933020089/OUTPUT/SeptemberFlorenceStreamflow_wDates.csv
- HYDRO_RST.2018-09-XX_00_00_DOMAIN1 -->
    test/input/florence_933020089/RESTART/ChannelRestart_2018-09-XX.csv

## Run the test
t-route simulation parameters for this test case are
stored as `.yaml` in `t-route/test/input/yaml/Florence_Benchmark.yaml`.

To execute the test:

```
bash

# navigate to python_routing_v02 directory
$ cd t-route/src/python_routing_v02

# compile the fortran Muskingum Cunge engine, the reservoir module,
# and the cythond code, then link all these together
# and import the result modules into the python site-packages
# for use in the script.
$ bash compiler.sh

# run the test with .yaml input
$ python compute_nhd_routing_SingleSeg_v02.py -f ../../test/input/yaml/Florence_Benchmark.yaml

```
## Test case data details
### `/DOMAIN`
This sub-directory contains essential model domain information,
specifically `LAKEPARM.nc` and `Route_Link.nc`, store
lake/reservoir and channel network parameters, respectively.
__Note:__ `LAKEPARM.nc` is not yet included.

### `/OUTPUT`
This sub-directory stores .csv files of WRF-Hydro simulated
lateral inflow and streamflow data.

Lateral inflow data are stored in `SeptemberFlorenceQlateral_wDates.csv`,
which contains the lateral inflow data with column
headers as timestamps.
Each row of data pertains to a
specific LinkID (channel node) in `DOMAIN/Route_Link.nc`.
Lateral inflow data have units of cubic meters per second.
WRF-simulated lateral inflows are the dynamic boundary
conditions of the t-route model(s).

Streamflow data are stored in `SeptemberFlorenceStreamflow_wDates.csv`.
Column headers denote simulation timestamps and each row pertains
to a specific LinkID (channel node) in  in `DOMAIN/Route_Link.nc`.
Streamflow inflow data have units of cubic meters per second.
WRF-simulated streamflows provide a baseline against which we
can validate the t-route Muskingum-Cunge model (i.e. t-roure
Muskingum-Cunge should generate the same results).

### `/RESTART`
This directory contains several `ChannelRestart*.csv`
files processed from the
WRF-Hydro `HYDRO_RST.*` files, which are used by
t-route to initialize flow and depth
states prior to a simulation. Streamflow data are in units of
cubic meters per second and depth data are in units of meters.
The `ChannelRestart*.csv` files contain only the flow
and depth data needed to initialize a model simulation. The contents
of `ChannelRestart*.csv` files are created by
`nhd_io.get_channel_restart_from_wrf_hydro`, which natively
reads-in the contents of a `HYDRO_RST.*` file.
