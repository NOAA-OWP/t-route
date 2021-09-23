#!/bin/bash -e


repo_dir=t-route/test/
testing_dir=/media/sf_Inland_hydrofabric/testing
huc=01a
nts=288

# get q_lat timeseries data

# run t-route
time python3 $repo_dir/refactor_to_t-route.py -huc $huc -nts $nts -d $testing_dir &>> $testing_dir/huc_"$huc".log

# plot metrics and find optimal parameters
python3 $srcDir/refactor_calibration.py -huc $huc -out $input_nwm_flows