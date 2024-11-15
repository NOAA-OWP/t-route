# Divergence Integration Test

### Author: Tadd Bindas

### Reason for this test:
There is a defined requirement for T-Route to be able to divert flows on river flowline breakages at man-made structures. One of these structures is the Old River Lock Control in Louisiana on the Mississippi River and Atchafalaya River confluence. This structure has USGS observed flow, but since there is no T-Route connectivity the flow is not taken into account

### What is tested:
The following config setting is evaluated:
```yaml
compute_parameters:
    data_assimilation_parameters:
        divergence_outflow : divergence_obs/old_river_lock.csv
```
which will remove observed flow from an inflow point (currently hard-coded), and drop the flow at an outflow point (beyond the water transer). 

The python functions involved are inside the `DataAssimilation.py` file and the `Divegence` class.

### How to run the test
You can either step through the provided Jupyter Notebook, which gives a detailed explaination about experiment set up and outputs, or run the following command from inside this directory:
```shell
python -m nwm_routing -V4 -f config.yaml
```
