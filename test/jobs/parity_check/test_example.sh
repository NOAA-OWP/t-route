#Remove outputs from previous runs.
rm /glade/u/home/jhreha/t-route/test/jobs/parity_check/chanobs_output/example1
rm /glade/u/home/jhreha/t-route/test/jobs/parity_check/chanobs_output/example2

#Each yaml file needs to have test outputs turned on and have a destination folder.
# Example
#output_parameters:
#    test_output: ../../t-route/test/jobs/parity_check/chanobs_output/example1
#    chanobs_output: 
#         chanobs_output_directory: ../../t-route/test/jobs/parity_check/chanobs_output/chanobs_output_test.ncdf

#Creates master branch outputs to use for benchmark.
# git checkout master #master
python -m nwm_routing -f -V3 parity_florence_example1.yaml 

#Creates test branch outputs for comparison against benchmark. 
# git checkout adam/make_chanobs #test branch
python -m nwm_routing -f -V3 parity_florence_example2.yaml 

#Saves a bash variable of the output directories for each file.
yaml1="$(cat parity_florence_example1.yaml | shyaml get-value output_parameters.test_output)"
yaml2="$(cat parity_florence_example2.yaml | shyaml get-value output_parameters.test_output)"

# #Runs the hydrotools Nash–Sutcliffe model efficiency coefficient tests to check for acceptable parity levels. 
python parity_NSE_example.py $yaml2 $yaml1