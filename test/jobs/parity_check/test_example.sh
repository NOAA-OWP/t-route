#Specify yaml1 file, will be used to generate the output of the master repository.
#Specify yaml2 file, will be used to generate the output of the test branch. 
yaml1='parity_florence_example1.yaml'
yaml2='parity_florence_example2.yaml'

#Saves a bash variable of the output directories for each file.
#Will be used to delete the output from the previous run and as an input argument for the parity check. 
yaml1_output="$(cat $yaml1 | shyaml get-value output_parameters.test_output)"
yaml2_output="$(cat $yaml2 | shyaml get-value output_parameters.test_output)"
#Remove outputs from previous runs.
rm $yaml1_output
rm $yaml2_output

#Each yaml file needs to have test outputs turned on and have a destination folder.
# Example
#output_parameters:
#    test_output: ../../t-route/test/jobs/parity_check/chanobs_output/example1
#    chanobs_output: 
#         chanobs_output_directory: ../../t-route/test/jobs/parity_check/chanobs_output/chanobs_output_test.ncdf

#Creates master branch outputs to use for benchmark.
git checkout master #master
python -m nwm_routing -f -V3 $yaml1 

#Creates test branch outputs for comparison against benchmark. 
git checkout adam/make_chanobs #test branch
# Please specify a test yaml used to run the test branch code. 
python -m nwm_routing -f -V3 $yaml2


# #Runs the hydrotools Nash–Sutcliffe model efficiency coefficient tests to check for acceptable parity levels. 
python parity_NSE_example.py $yaml1_output $yaml2_output