#Remove outputs from previous runs.
rm /glade/u/home/jhreha/t-route/test/jobs/parity_check/chanobs_output/example1
rm /glade/u/home/jhreha/t-route/test/jobs/parity_check/chanobs_output/example2
#Creates master branch outputs to use for benchmark.
git checkout master #master
python -m nwm_routing -f -V3 parity_florence_example1.yaml 
#Creates test branch outputs for comparison against benchmark. 
git checkout adam/make_chanobs #test branch
python -m nwm_routing -f -V3 parity_florence_example2.yaml 
#Runs the hydrotools Nash–Sutcliffe model efficiency coefficient tests to check for acceptable parity levels. 
python parity_NSE_generalized.py