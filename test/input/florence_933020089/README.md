# Hurricane Florence Test Case
This directory contains data from a WRF-Hydro simulation of Hurricane Florence on a domain in North Carolina, USA. Because WRF-Hydro output data can be quite bulky, the lateral inflow and streamflow data found here have been constructed by pre-processing contents of many .CHRTOUT files. 

## Run the test
t-route simulation parameters for this test case are stored as `.yaml` in `t-route/test/input/yaml/Florence_Benchmark.yaml`. To run execite the test:

```bash
# navigate to python_routing_v02 directory
$ cd t-route/src/python_routing_v02 

# run the test with .yaml input
$ python compute_nhd_routing_SingleSeg_v02.py -f ../../test/input/yaml/Florence_Benchmark.yaml

```
## Test case data details
### `/DOMAIN`
This sub-directory contains essential model domain information, specifically `LAKEPARM.nc` and `Route_Link.nc`, store lake/reservoir and channel network parameters, respectively. 

### `/OUTPUT`
This sub-directory stores .csv files of WRF-Hydro simulated lateral inflow and streamflow data. 

Lateral inflow data are stored in `SeptemberFlorenceQlateral_wDates.csv` and `SeptemberFlorenceQlateral.csv`. Each file contains the lateral inflow data, however column headers are timestamps in `SeptemberFlorenceQlateral_wDates.csv` and simply timestep integers (0, 1, ..., n) in `SeptemberFlorenceQlateral.csv`. Each row of data pertains to a specific LinkID (channel node) in `DOMAIN/Route_Link.nc`. Lateral inflow data have units of cubic meters per second. WRF-simulated lateral inflows are the dynamiuc boundary conditions of the t-route model(s).

Streamflow data are stored in `SeptemberFlorenceStreamflow_wDates.csv`. Column headers denote simulation timestamps and each row pertains to a specific LinkID (channel node) in  in `DOMAIN/Route_Link.nc`. Streamflow inflow data have units of cubic meters per second. WRF-simulated streamflows provide a baseline against which we can validate the t-route Muskingum-Cunge model (i.e. t-roure Muskingum-Cunge should generate the same results).

### `/RESTART`
This directory contains several WRF-Hydro `HYDRO_RST.*` files, which are used by t-route to initialize flow and depth states prior to a simulation. Streamflow data are in units of cubic meters per second and depth data are in units of meters. Additionally, `ChannelRestart*.csv` files contain only the flow and depth data needed to initialize a model simulation. The contents of `ChannelRestart*.csv` files are created by `nhd_io.get_stream_restart_from_wrf_hydro`, which reads-in the contents of a `HYDRO_RST.*` file.