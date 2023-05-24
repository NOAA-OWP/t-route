

import sys
import numpy as np
import pandas as pd

sys.path.append("src/")
import bmi_troute
import bmi_reservoirs

#----------------------------------------------------------------
# Full model
#----------------------------------------------------------------
# Initialize model with configuration file
full_model = bmi_troute.bmi_troute()
full_model.initialize(bmi_cfg_file='test/BMI/bmi_example.yaml')

# Set static values
#Segment parameters
full_model.set_value('segment_id', np.array([1056,385,156,158,159,157]))
full_model.set_value('segment_toid', np.array([157,157,157,159,160,158]))
full_model.set_value('dx', np.array([3383.185374,7247.358765,2090.045820,2059.280869,3909.862245,1568.900922]))
full_model.set_value('n', np.array([0.060000,0.059635,0.050000,0.057287,0.050000,0.055400]))
full_model.set_value('ncc', np.array([0.120000,0.119270,0.100000,0.114575,0.100000,0.110800]))
full_model.set_value('s0', np.array([0.011690,0.018579,0.000211,0.030771,0.002016,0.017703]))
full_model.set_value('bw', np.array([2.860991,2.809054,19.379018,6.677849,20.134363,9.857479]))
full_model.set_value('tw', np.array([4.768318,4.681756,32.298364,11.129749,33.557271,16.429131]))
full_model.set_value('twcc', np.array([14.304953,14.045269,96.895086,33.389247,100.671812,49.287396]))
full_model.set_value('alt', np.array([1.0,1.0,1.0,1.0,1.0,1.0]))
full_model.set_value('musk', np.array([3600.0,3600.0,3600.0,3600.0,3600.0,3600.0]))
full_model.set_value('musx', np.array([0.2,0.2,0.2,0.2,0.2,0.2]))
full_model.set_value('cs', np.array([0.585855,0.610357,0.251915,0.581235,0.247701,0.519852]))

#Waterbody parameters
full_model.set_value('waterbody_id', np.array([157]))
full_model.set_value('waterbody_toid', np.array([158]))
full_model.set_value('LkArea', np.array([61.150299]))
full_model.set_value('LkMxE', np.array([201.179993]))
full_model.set_value('OrificeA', np.array([1.0]))
full_model.set_value('OrificeC', np.array([0.1]))
full_model.set_value('OrificeE', np.array([190.559998]))
full_model.set_value('WeirC', np.array([0.4]))
full_model.set_value('WeirE', np.array([199.586993]))
full_model.set_value('WeirL', np.array([10.0]))
full_model.set_value('ifd', np.array([0.9]))
full_model.set_value('reservoir_type', np.array([2]))

#-----------------------
# Split model
#-----------------------
# Initialize routing models with configuration file
upper_routing_model = bmi_troute.bmi_troute()
upper_routing_model.initialize(bmi_cfg_file='test/BMI/bmi_upper_example.yaml')

lower_routing_model = bmi_troute.bmi_troute()
lower_routing_model.initialize(bmi_cfg_file='test/BMI/bmi_lower_example.yaml')

# Initialize reservoir model
reservoir_model = bmi_reservoirs.bmi_reservoir()
reservoir_model.initialize(bmi_cfg_file='test/BMI/bmi_reservoir_example.yaml')

# Set static values
#Segment parameters, upper
upper_routing_model.set_value('segment_id', np.array([1056,385,156]))
upper_routing_model.set_value('segment_toid', np.array([157,157,157]))
upper_routing_model.set_value('dx', np.array([3383.185374,7247.358765,2090.045820]))
upper_routing_model.set_value('n', np.array([0.060000,0.059635,0.050000]))
upper_routing_model.set_value('ncc', np.array([0.120000,0.119270,0.100000]))
upper_routing_model.set_value('s0', np.array([0.011690,0.018579,0.000211]))
upper_routing_model.set_value('bw', np.array([2.860991,2.809054,19.379018]))
upper_routing_model.set_value('tw', np.array([4.768318,4.681756,32.298364]))
upper_routing_model.set_value('twcc', np.array([14.304953,14.045269,96.895086]))
upper_routing_model.set_value('alt', np.array([1.0,1.0,1.0]))
upper_routing_model.set_value('musk', np.array([3600.0,3600.0,3600.0]))
upper_routing_model.set_value('musx', np.array([0.2,0.2,0.2]))
upper_routing_model.set_value('cs', np.array([0.585855,0.610357,0.251915]))

#Segment parameters, lower
lower_routing_model.set_value('segment_id', np.array([158,159,157]))
lower_routing_model.set_value('segment_toid', np.array([159,160,158]))
lower_routing_model.set_value('dx', np.array([2059.280869,3909.862245,1568.900922]))
lower_routing_model.set_value('n', np.array([0.057287,0.050000,0.055400]))
lower_routing_model.set_value('ncc', np.array([0.114575,0.100000,0.110800]))
lower_routing_model.set_value('s0', np.array([0.030771,0.002016,0.017703]))
lower_routing_model.set_value('bw', np.array([6.677849,20.134363,9.857479]))
lower_routing_model.set_value('tw', np.array([11.129749,33.557271,16.429131]))
lower_routing_model.set_value('twcc', np.array([33.389247,100.671812,49.287396]))
lower_routing_model.set_value('alt', np.array([1.0,1.0,1.0]))
lower_routing_model.set_value('musk', np.array([3600.0,3600.0,3600.0]))
lower_routing_model.set_value('musx', np.array([0.2,0.2,0.2]))
lower_routing_model.set_value('cs', np.array([0.581235,0.247701,0.519852]))

#Waterbody parameters
reservoir_model.set_value('waterbody_id', np.array([157]))
#reservoir_model.set_value('waterbody_toid', np.array([158]))
reservoir_model.set_value('LkArea', np.array([61.150299]))
reservoir_model.set_value('LkMxE', np.array([201.179993]))
reservoir_model.set_value('OrificeA', np.array([1.0]))
reservoir_model.set_value('OrificeC', np.array([0.1]))
reservoir_model.set_value('OrificeE', np.array([190.559998]))
reservoir_model.set_value('WeirC', np.array([0.4]))
reservoir_model.set_value('WeirE', np.array([199.586993]))
reservoir_model.set_value('WeirL', np.array([10.0]))
reservoir_model.set_value('ifd', np.array([0.9]))
reservoir_model.set_value('reservoir_type', np.array([2]))
reservoir_model.set_value('lake_surface__elevation', np.array([-1.000000e+09]))

# Create random forcing values
forcing_vals = np.random.gamma(1,1,6*120).reshape(120,6)
gage_vals = np.random.gamma(1,1,1*120)
gage_times = np.arange(0.0,120*3600,3600)
# input some NaNs to see how persistence module handles missing values
gage_vals[25:73] = np.nan
reservoir_model.set_value('gage_observations', gage_vals)
reservoir_model.set_value('gage_time', gage_times)

full_fvd = pd.DataFrame()
full_wbdy = pd.DataFrame()
upper_fvd = pd.DataFrame()
lower_fvd = pd.DataFrame()
res_fvd = pd.DataFrame()

update_until_duration = 3600 # run for 1 hour
for hr in range(120):
    # Full network
    # Set dynamic values
    full_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,:])
    #full_model.set_value('gage_observations', gage_vals)
    #full_model.set_value('gage_time', gage_times)

    full_model.update_until(update_until_duration)

    flowveldepth = pd.DataFrame(full_model.get_value('fvd_results').reshape(6,36),index=full_model.get_value('fvd_index'))
    full_fvd = pd.concat([full_fvd, flowveldepth], axis=1)
    
    # Split network
    # Set dynamic values
    upper_routing_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,0:3])
    
    upper_routing_model.update_until(update_until_duration)

    fvd_index = upper_routing_model.get_value('fvd_index')
    fvd_num_columns = int(update_until_duration/upper_routing_model.get_time_step()*3) # duration of update divided by model time steps times 3 variables (fvd)
    
    flowveldepth = pd.DataFrame(upper_routing_model.get_value('fvd_results').reshape(len(fvd_index),fvd_num_columns),index=fvd_index)
    upper_fvd = pd.concat([upper_fvd, flowveldepth], axis=1)
    
    reservoir_inflow = flowveldepth.loc[[1056,385,156]].sum()[0::3].to_numpy()
    reservoir_model.set_value('lake_water~incoming__volume_flow_rate', reservoir_inflow)

    reservoir_model.update_until(update_until_duration)
    
    offnetwork_upstream_fvd = np.asarray(
        [
        item for sublist in zip(reservoir_model._model._outflow_list,
                                np.zeros(len(reservoir_model._model._outflow_list)).tolist(),
                                reservoir_model._model._water_elevation_list) for item in sublist
        ]
    )
    reservoir_model._model._outflow_list = []
    reservoir_model._model._water_elevation_list = []
    
    lower_routing_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,3:6])
    lower_routing_model.set_value('upstream_id', np.array([int(157)]))
    lower_routing_model.set_value('upstream_fvd', offnetwork_upstream_fvd)
    lower_routing_model.update_until(3600)

    fvd_index = lower_routing_model.get_value('fvd_index')
    fvd_num_columns = int(update_until_duration/lower_routing_model.get_time_step()*3) # duration of update divided by model time steps times 3 variables (fvd)
    
    flowveldepth = pd.DataFrame(lower_routing_model.get_value('fvd_results').reshape(len(fvd_index),fvd_num_columns),index=fvd_index)    
    lower_fvd = pd.concat([lower_fvd, flowveldepth], axis=1)
    res_fvd = pd.concat([res_fvd, pd.DataFrame(offnetwork_upstream_fvd.reshape(12,3))], axis=0)
