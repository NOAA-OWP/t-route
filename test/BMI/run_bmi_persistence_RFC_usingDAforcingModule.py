#import sys
import numpy as np
import pandas as pd

#sys.path.append("src/")
import bmi_troute
import bmi_reservoirs
import bmi_DAforcing

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
full_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_large_example.yaml')

# Set static values
#flowpaths layer parameters
full_model.set_value('id', np.array([1056,385,156,158,159,157]))
full_model.set_value('toid', np.array([157,157,157,159,160,158]))
full_model.set_value('lengthkm', np.array([3.383,7.247,2.090,2.059,3.909,1.568]))
#flowpath_attributes layer parameters
full_model.set_value('attributes_id', np.array([1056,385,156,158,159,157]))
full_model.set_value('rl_gages', np.array(['08117995', '08124000', '08130700', 'NULL', 'NULL', 'NULL']))
full_model.set_value('rl_NHDWaterbodyComID', np.array(['NULL', 'NULL', 'NULL','NULL', 'NULL', 157 ]))
full_model.set_value('n', np.array([0.060000,0.059635,0.050000,0.057287,0.050000,0.055400]))
full_model.set_value('nCC', np.array([0.120000,0.119270,0.100000,0.114575,0.100000,0.110800]))
full_model.set_value('So', np.array([0.011690,0.018579,0.000211,0.030771,0.002016,0.017703]))
full_model.set_value('BtmWdth', np.array([2.860991,2.809054,19.379018,6.677849,20.134363,9.857479]))
full_model.set_value('TopWdth', np.array([4.768318,4.681756,32.298364,11.129749,33.557271,16.429131]))
full_model.set_value('TopWdthCC', np.array([14.304953,14.045269,96.895086,33.389247,100.671812,49.287396]))
full_model.set_value('alt', np.array([1.0,1.0,1.0,1.0,1.0,1.0]))
full_model.set_value('MusK', np.array([3600.0,3600.0,3600.0,3600.0,3600.0,3600.0]))
full_model.set_value('MusX', np.array([0.2,0.2,0.2,0.2,0.2,0.2]))
full_model.set_value('ChSlp', np.array([0.585855,0.610357,0.251915,0.581235,0.247701,0.519852]))
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
upper_routing_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_upper_example.yaml')

lower_routing_model = bmi_troute.bmi_troute()
lower_routing_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_lower_example.yaml')

# Initialize reservoir model
reservoir_model = bmi_reservoirs.bmi_reservoir()
reservoir_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_reservoir_example.yaml')

# Set static values
#flowpaths layer parameters
upper_routing_model.set_value('id', np.array(['wb-1056','wb-385','wb-156']))
upper_routing_model.set_value('toid', np.array(['nex-157','nex-157','nex-157']))
upper_routing_model.set_value('lengthkm', np.array([3.383,7.247,2.090]))
#flowpath_attributes layer parameters
upper_routing_model.set_value('attributes_id', np.array(['wb-1056','wb-385','wb-156']))
upper_routing_model.set_value('rl_gages', np.array(['08117995', '08124000', '08130700']))
upper_routing_model.set_value('rl_NHDWaterbodyComID', np.array([None, None, None]))
upper_routing_model.set_value('n', np.array([0.060000,0.059635,0.050000]))
upper_routing_model.set_value('nCC', np.array([0.120000,0.119270,0.100000]))
upper_routing_model.set_value('So', np.array([0.011690,0.018579,0.000211]))
upper_routing_model.set_value('BtmWdth', np.array([2.860991,2.809054,19.379018]))
upper_routing_model.set_value('TopWdth', np.array([4.768318,4.681756,32.298364]))
upper_routing_model.set_value('TopWdthCC', np.array([14.304953,14.045269,96.895086]))
upper_routing_model.set_value('alt', np.array([1.0,1.0,1.0]))
upper_routing_model.set_value('MusK', np.array([3600.0,3600.0,3600.0]))
upper_routing_model.set_value('MusX', np.array([0.2,0.2,0.2]))
upper_routing_model.set_value('ChSlp', np.array([0.585855,0.610357,0.251915]))
upper_routing_model.set_value('network_id', np.array(['wb-1056','wb-385','wb-156']))
upper_routing_model.set_value('hydroseq', np.array([1,2,3]) )
upper_routing_model.set_value('hl_uri', np.array(['Gages-08117995', 'Gages-08124000', 'Gages-08130700']))

#flowpaths layer parameters, lower
lower_routing_model.set_value('id', np.array(['wb-158','wb-159','wb-157']))
lower_routing_model.set_value('toid', np.array(['nex-159','nex-160','nex-158']))
lower_routing_model.set_value('lengthkm', np.array([2.059,3.909,1.568]))
#flowpath_attributes parameters
lower_routing_model.set_value('attributes_id', np.array(['wb-158','wb-159','wb-157']))
lower_routing_model.set_value('rl_gages', np.array([None, None, None]))
lower_routing_model.set_value('rl_NHDWaterbodyComID', np.array([None, None, None]))
lower_routing_model.set_value('n', np.array([0.057287,0.050000,0.055400]))
lower_routing_model.set_value('nCC', np.array([0.114575,0.100000,0.110800]))
lower_routing_model.set_value('So', np.array([0.030771,0.002016,0.017703]))
lower_routing_model.set_value('BtmWdth', np.array([6.677849,20.134363,9.857479]))
lower_routing_model.set_value('TopWdth', np.array([11.129749,33.557271,16.429131]))
lower_routing_model.set_value('TopWdthCC', np.array([33.389247,100.671812,49.287396]))
lower_routing_model.set_value('alt', np.array([1.0,1.0,1.0]))
lower_routing_model.set_value('MusK', np.array([3600.0,3600.0,3600.0]))
lower_routing_model.set_value('MusX', np.array([0.2,0.2,0.2]))
lower_routing_model.set_value('ChSlp', np.array([0.581235,0.247701,0.519852]))
lower_routing_model.set_value('network_id', np.array(['wb-158','wb-159','wb-157']))
lower_routing_model.set_value('hydroseq', np.array([5, 6, 4]) )
lower_routing_model.set_value('hl_uri', np.array([None, None, None]))


# build RFC DA forcing and related inputs
DAforcing_model = bmi_DAforcing.bmi_DAforcing()
DAforcing_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_large_example.yaml')
# input argument to DA forcing preprocess module
#DAforcing_model.set_value('rfc_gage_id', np.array(['KNFC1']))
DAforcing_model.set_value('waterbody__lake_number', np.array(['347987'])) # assume 157 <- 347987
DAforcing_model.set_value('waterbody__type_number', np.array([2]))
DAforcing_model.update()


# Initialize reservoir model
reservoir_model = bmi_reservoirs.bmi_reservoir()
reservoir_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_reservoir_example.yaml')

# move DA forcing data obtained from DAforcing_model to reservoir_model
#Waterbody parameters
reservoir_model.set_value('hl_link', np.array([347987]))
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
reservoir_model.set_value('reservoir_type', np.array(DAforcing_model._model._res_type))

if DAforcing_model._model._res_type==4 or DAforcing_model._model._res_type==5:
    reservoir_model.set_value('lake_number', DAforcing_model._values['lake_number'] )
    reservoir_model.set_value('reservoir_type', np.array(DAforcing_model._model._res_type))
    reservoir_model.set_value('use_RFC', np.array(DAforcing_model._model._use_RFC)) 
    reservoir_model.set_value('rfc_timeseries_discharges',  np.array(DAforcing_model._model._rfc_timeseries_discharges))
    reservoir_model.set_value('rfc_timeseries_idx',  np.array(DAforcing_model._model._rfc_timeseries_idx)) 
    reservoir_model.set_value('rfc_timeseries_update_time',  np.array(DAforcing_model._model._rfc_update_time))
    reservoir_model.set_value('rfc_da_time_step',  np.array(DAforcing_model._model._da_time_step))
    reservoir_model.set_value('rfc_total_counts',  np.array(DAforcing_model._model._total_counts))  
    reservoir_model.set_value('rfc_timeseries_file',  np.array(DAforcing_model._model._rfc_timeseries_file))
    reservoir_model.set_value('lake_surface__elevation', np.array([-1.000000e+09]))

if DAforcing_model._model._res_type==2 or DAforcing_model._model._res_type==3:
    #DA forcing for usgs_df
    upper_routing_model.set_value('gages', np.array(['1056','08117995','385','08124000','156','08130700']))
    #upper_routing_model.set_value('usgs_timeslice_discharge', np.array([1,2,3,4,5,6]))
    #upper_routing_model.set_value('usgs_timeslice_stationId', np.array(['01012570', '01013500', '01012515', '01012570', '01013500', '01012515']))
    #upper_routing_model.set_value('usgs_timeslice_time', np.array(['2015-08-16_00:00:00','2015-08-16_00:03:00','2015-08-16_00:03:00',
    #                                                '2015-08-16_00:15:00', '2015-08-16_00:18:00','2015-08-16_00:18:00']))
    #upper_routing_model.set_value('usgs_timeslice_discharge_quality', np.array([100,100,100,100,100,100]))
    upper_routing_model.set_value('usgs_timeslice_discharge', np.array(DAforcing_model._model._timeslice_obs))
    upper_routing_model.set_value('usgs_timeslice_stationId', np.array(DAforcing_model._model._timeslice_stationId))
    upper_routing_model.set_value('usgs_timeslice_time', np.array(DAforcing_model._model._timeslice_datetime))
    upper_routing_model.set_value('usgs_timeslice_discharge_quality', np.array(DAforcing_model._model._timeslice_discharge_quality))

    #DA forcing for last_obs_df
    upper_routing_model.set_value('lastobs_discharge', np.array(DAforcing_model._model._lastobs_discharge))
    upper_routing_model.set_value('time_since_lastobs', np.array(DAforcing_model._model._time_since_lastobs))
    upper_routing_model.set_value('lastobs_stationId', np.array(DAforcing_model._model._lastobs_stationId))    

# Create random forcing values
forcing_vals = np.random.gamma(1,1,6*120).reshape(120,6)
gage_vals = np.random.gamma(1,1,1*120)
gage_times = np.arange(0.0,120*3600,3600)
# input some NaNs to see how persistence module handles missing values
gage_vals[25:73] = np.nan

reservoir_model.set_value('gage_observations', gage_vals)
reservoir_model.set_value('gage_time', gage_times)
# input for RFC reservoir DA

#reservoir_model.set_value('rfc_timeseries_offset_hours', np.array([28]))
#reservoir_model.set_value('rfc_forecast_persist_days',np.array([11]))
        
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
    /home/dongha.kim/github/t-route/test/inputupper_routing_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,0:3])
    upper_routing_model.update_until(3600)
    flowveldepth, wbdy = create_output_dataframes(upper_routing_model._model._run_results, 
                                                  3600/upper_routing_model._model._time_step, 
                                                  upper_routing_model._model._network.waterbody_dataframe)
    upper_fvd = pd.concat([upper_fvd, flowveldepth], axis=1)

    reservoir_inflow = flowveldepth.loc[[1056,385,156]].sum()[0::3].to_numpy()
    # test
    '''
    reservoir_inflow = np.array([0.038284,
                                 0.105146,
                                 0.219112,
                                 0.381832,
                                 0.544385,
                                 0.701310,
                                 0.849002,
                                 0.985553,
                                 1.111356,
                                 1.224970,
                                 1.326913,
                                 1.419024])
    '''
    reservoir_inflow = np.array([187,
                                185,
                                183,
                                181,
                                179,
                                177,
                                175,
                                173,
                                171,
                                169,
                                167,
                                165,
                                163,
                                161,
                                159,
                                157,
                                155,
                                153,
                                151,
                                149,
                                147,
                                145,
                                143,
                                141,
                                139])

    reservoir_model.set_value('lake_water~incoming__volume_flow_rate', reservoir_inflow)
    reservoir_model.update_until(3600)
    ## here here lover
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

full_fvd.to_csv('/home/dongha.kim/github/t-route/temp_output/persistence/full_fvd.csv')
upper_fvd.to_csv('/home/dongha.kim/github/t-route/temp_output/persistence/upper_fvd.csv')
lower_fvd.to_csv('/home/dongha.kim/github/t-route/temp_output/persistence/lower_fvd.csv')
res_fvd.to_csv('/home/dongha.kim/github/t-route/temp_output/persistence/res_fvd.csv')
pd.DataFrame({'Gage_Time': gage_times,
              'Gage_Values': gage_vals}).to_csv('/home/dongha.kim/github/t-route/temp_output/persistence/gage.csv')
