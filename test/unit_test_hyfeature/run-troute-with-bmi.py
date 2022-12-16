import sys
sys.path.append("../../src/")
import bmi_troute # This is the BMI t-route that we will be running from the file: bmi_troute.py 
import troute.main_utilities as mu # This is used to read q_lateral data from files, won't be needed when t-route bmi is run with model engine
import troute.nhd_preprocess as nhd_prep # This is used to read q_lateral data from files, won't be needed when t-route bmi is run with model engine
model = bmi_troute.bmi_troute()

# this call covers lines 41 ~ 110 of def main_v04
model.initialize(bmi_cfg_file='unittest_hyfeature.yaml')

# this call overs lines 112 ~ 126
(run_sets, da_sets, parity_sets) = mu.build_run_sets(model._network,
                                                     model._supernetwork_parameters,
                                                     model._forcing_parameters,
                                                     model._compute_parameters, 
                                                     model._data_assimilation_parameters,
                                                     model._output_parameters,
                                                     model._parity_parameters)

# retrieve forcing input data such as lateral flow and coastal boundary data
cpu_pool = model._compute_parameters.get("cpu_pool", None)
model._network.assemble_forcings(run_sets[0], model._forcing_parameters, model._hybrid_parameters, cpu_pool)
import pdb; pdb.set_trace()
# move lateral flow data into BMI compatible variable (model._values)
model.set_value('land_surface_water_source__volume_flow_rate', model._network._qlateral.to_numpy())
model.set_value('coastal_boundary__depth', model._network._coastal_boundary_depth_df.to_numpy())

i=2
#model.update()