import sys
sys.path.append("src/")
import bmi_troute # This is the BMI t-route that we will be running from the file: bmi_troute.py 
import troute.main_utilities as mu # This is used to read q_lateral data from files, won't be needed when t-route bmi is run with model engine
import troute.nhd_preprocess as nhd_prep # This is used to read q_lateral data from files, won't be needed when t-route bmi is run with model engine
model = bmi_troute.bmi_troute()
import pdb; pdb.set_trace()
#model.initialize(bmi_cfg_file='test/BMI/test_AnA_bmi.yaml')
model.initialize(bmi_cfg_file='test/unit_test_hyfeature/unittest_hyfeature.yaml')
import pdb; pdb.set_trace()
(run_sets, da_sets, parity_sets) = mu.build_run_sets(model._network,
                                                     model._forcing_parameters,
                                                     model._compute_parameters, 
                                                     model._data_assimilation_parameters,
                                                     model._output_parameters,
                                                     model._parity_parameters)


q_lateral, coastal_boundary_depth_df = nhd_prep.nhd_forcing(run_sets[0], 
                                                            model._forcing_parameters, 
                                                            model._hybrid_parameters,
                                                            model._network.segment_index,
                                                            model._compute_parameters.get('cpu_pool', None),
                                                            model._network._t0,             
                                                            model._network._coastal_boundary_depth_df,
        )
import pdb; pdb.set_trace()

model.set_value('land_surface_water_source__volume_flow_rate', q_lateral.to_numpy())


