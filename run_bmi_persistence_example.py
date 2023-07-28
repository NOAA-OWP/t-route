
#import sys
import numpy as np
import pandas as pd

#sys.path.append("src/")
import bmi_troute
import bmi_reservoirs

def create_output_dataframes(results, nts, waterbodies_df):
    """
    Run this model into the future, updating the state stored in the provided model dict appropriately.
    Note that the model assumes the current values set for input variables are appropriately for the time
    duration of this update (i.e., ``dt``) and do not need to be interpolated any here.
    Parameters
    ----------
    results: list
        The results from nwm_routing.
    nts: int
        The number of time steps the model was run.
    waterbodies_df: pd.DataFrame
        Dataframe containing waterbody parameters (specifically, IDs stored in index)
    link_lake_crosswalk: dict #TODO: Can we remove this?
        Relates lake ids to outlet link ids.
    Returns
    -------
    q_channel_df: pandas.core.series.Series
        Streamflow rate for each segment
    v_channel_df: pandas.core.series.Series
        Streamflow velocity for each segment
    d_channel_df: pandas.core.series.Series
        Streamflow depth for each segment
    i_lakeout_df: pandas.core.series.Series
        Inflow for each waterbody
    q_lakeout_df: pandas.core.series.Series
        Outflow for each waterbody
    d_lakeout_df: pandas.core.series.Series
        Water elevation for each waterbody
    """
    qvd_columns = pd.MultiIndex.from_product(
        [range(int(nts)), ["q", "v", "d"]]
    ).to_flat_index()
    
    flowveldepth = pd.concat(
        [pd.DataFrame(r[1], index=r[0], columns=qvd_columns) for r in results], copy=False,
    )
    
    # create waterbody dataframe for output to netcdf file
    i_columns = pd.MultiIndex.from_product(
        [range(int(nts)), ["i"]]
    ).to_flat_index()
    
    wbdy = pd.concat(
        [pd.DataFrame(r[6], index=r[0], columns=i_columns) for r in results],
        copy=False,
    )
    
    wbdy_id_list = waterbodies_df.index.values.tolist()

    i_lakeout_df = wbdy.loc[wbdy_id_list].iloc[:,-1]
    q_lakeout_df = flowveldepth.loc[wbdy_id_list].iloc[:,-3]
    d_lakeout_df = flowveldepth.loc[wbdy_id_list].iloc[:,-1]
    #lakeout = pd.concat([i_df, q_df, d_df], axis=1)
    
    # replace waterbody lake_ids with outlet link ids
    #TODO Update the following line to fit with HyFeatures. Do we need to replace IDs? Or replace
    # waterbody_ids with the downstream segment?
    #flowveldepth = _reindex_lake_to_link_id(flowveldepth, link_lake_crosswalk)
    
    q_channel_df = flowveldepth.iloc[:,-3]
    v_channel_df = flowveldepth.iloc[:,-2]
    d_channel_df = flowveldepth.iloc[:,-1]
    
    segment_ids = flowveldepth.index.values.tolist()

    return flowveldepth, wbdy #q_channel_df, v_channel_df, d_channel_df, i_lakeout_df, q_lakeout_df, d_lakeout_df#, wbdy_id_list, 

#----------------------------------------------------------------
# Full model
#----------------------------------------------------------------
# Initialize model with configuration file
full_model = bmi_troute.bmi_troute()
full_model.initialize(bmi_cfg_file='../test/BMI/bmi_example.yaml')

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
reservoir_model.initialize()

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
for hr in range(120):
    '''
    # Full network
    # Set dynamic values
    full_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,:])
    full_model.set_value('gage_observations', gage_vals)
    full_model.set_value('gage_time', gage_times)
    import pdb; pdb.set_trace()
    full_model.update_until(3600)
    
    flowveldepth, wbdy = create_output_dataframes(full_model._model._run_results, 
                                                  3600/full_model._model._time_step, 
                                                  full_model._model._network.waterbody_dataframe)
    pdb.set_trace()
    full_fvd = pd.concat([full_fvd, flowveldepth], axis=1)
    full_wbdy = pd.concat([full_wbdy, wbdy], axis=1)
    '''
    # Split network
    # Set dynamic values
    upper_routing_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,0:3])

    upper_routing_model.update_until(3600)

    flowveldepth, wbdy = create_output_dataframes(upper_routing_model._model._run_results, 
                                                  3600/upper_routing_model._model._time_step, 
                                                  upper_routing_model._model._network.waterbody_dataframe)
    upper_fvd = pd.concat([upper_fvd, flowveldepth], axis=1)
    
    reservoir_inflow = flowveldepth.loc[[1056,385,156]].sum()[0::3].to_numpy()
    reservoir_model.set_value('lake_water~incoming__volume_flow_rate', reservoir_inflow)

    reservoir_model.update_until(3600)
    
    upstream_fvd = np.asarray(
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
    lower_routing_model.set_value('upstream_fvd', upstream_fvd)
    lower_routing_model.update_until(3600)

    flowveldepth, wbdy = create_output_dataframes(lower_routing_model._model._run_results, 
                                                  3600/lower_routing_model._model._time_step, 
                                                  lower_routing_model._model._network.waterbody_dataframe)
    
    lower_fvd = pd.concat([lower_fvd, flowveldepth], axis=1)
    res_fvd = pd.concat([res_fvd, pd.DataFrame(upstream_fvd.reshape(12,3))], axis=0)
    
import pdb; pdb.set_trace()
full_fvd.to_csv('/home/sean.horvath/projects/data/test_reservoir_module/persistence/full_fvd.csv')
upper_fvd.to_csv('/home/sean.horvath/projects/data/test_reservoir_module/persistence/upper_fvd.csv')
lower_fvd.to_csv('/home/sean.horvath/projects/data/test_reservoir_module/persistence/lower_fvd.csv')
res_fvd.to_csv('/home/sean.horvath/projects/data/test_reservoir_module/persistence/res_fvd.csv')
pd.DataFrame({'Gage_Time': gage_times,
              'Gage_Values': gage_vals}).to_csv('/home/sean.horvath/projects/data/test_reservoir_module/persistence/gage.csv')

