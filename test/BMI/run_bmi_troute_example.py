#import sys
import numpy as np
import pandas as pd
import pickle

#sys.path.append("src/")
import bmi_troute
import bmi_DAforcing


def numeric_id(df):
    toid = df['toid'].split('-')[-1]
    df['toid'] = int(toid)
    return df

def read_bmi_input_data(filename):
    with open(filename, 'rb') as h:
        df, wbdy_df = pickle.load(h)

    #Remove string prefix from wbdy['toid']
    wbdy_df =pd.DataFrame(wbdy_df)
    wbdy_df = wbdy_df.apply(numeric_id, axis=1)

    df_dict = df.reset_index().to_dict('list')
    wbdy_dict = wbdy_df.reset_index().to_dict('list')

    return df_dict, wbdy_dict


###############################################
#TEMP
from datetime import timedelta
def lastobs_df_output(
    lastobs_df,
    dt,
    nts,
    t0,
    gages,
):

    # join gageIDs to lastobs_df
    lastobs_df = lastobs_df.join(gages)

    # timestamp of last simulation timestep
    modelTimeAtOutput = t0 + timedelta(seconds = nts * dt)
    modelTimeAtOutput_str = modelTimeAtOutput.strftime('%Y-%m-%d_%H:%M:%S')

    # timestamp of last observation
    var = [timedelta(seconds=d) for d in lastobs_df.time_since_lastobs.fillna(0)]
    lastobs_timestamp = [modelTimeAtOutput - d for d in var]
    lastobs_timestamp_str = [d.strftime('%Y-%m-%d_%H:%M:%S') for d in lastobs_timestamp]
    lastobs_timestamp_str_array = np.asarray(lastobs_timestamp_str,dtype = '|S19').reshape(len(lastobs_timestamp_str),1)

    return lastobs_timestamp_str_array
###################################################



#----------------------------------------------------------------
# Full model
#----------------------------------------------------------------
# Initialize model with configuration file
full_model = bmi_troute.bmi_troute()
full_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_large_example.yaml')

# Set static values
# Read from pickle file to imitate BMI set_value function
df_dict, wbdy_dict = read_bmi_input_data('/home/dongha.kim/github/t-route/temp_input/bmi_input/bmi_input_data.pickle')

#Segment parameters
full_model.set_value('segment_id', np.array(df_dict['key']))
full_model.set_value('segment_toid', np.array(df_dict['downstream']))
full_model.set_value('dx', np.array(df_dict['dx']))
full_model.set_value('n', np.array(df_dict['n']))
full_model.set_value('ncc', np.array(df_dict['ncc']))
full_model.set_value('s0', np.array(df_dict['s0']))
full_model.set_value('bw', np.array(df_dict['bw']))
full_model.set_value('tw', np.array(df_dict['tw']))
full_model.set_value('twcc', np.array(df_dict['twcc']))
full_model.set_value('alt', np.array(df_dict['alt']))
full_model.set_value('musk', np.array(df_dict['musk']))
full_model.set_value('musx', np.array(df_dict['musx']))
full_model.set_value('cs', np.array(df_dict['cs']))

#Waterbody parameters
full_model.set_value('waterbody_id', np.array(wbdy_dict['wb-id']))
full_model.set_value('waterbody_toid', np.array(wbdy_dict['toid']))
full_model.set_value('LkArea', np.array(wbdy_dict['LkArea']))
full_model.set_value('LkMxE', np.array(wbdy_dict['LkMxE']))
full_model.set_value('OrificeA', np.array(wbdy_dict['OrificeA']))
full_model.set_value('OrificeC', np.array(wbdy_dict['OrificeC']))
full_model.set_value('OrificeE', np.array(wbdy_dict['OrificeE']))
full_model.set_value('WeirC', np.array(wbdy_dict['WeirC']))
full_model.set_value('WeirE', np.array(wbdy_dict['WeirE']))
full_model.set_value('WeirL', np.array(wbdy_dict['WeirL']))
full_model.set_value('ifd', np.array(wbdy_dict['ifd']))
full_model.set_value('reservoir_type', np.array([1,1,1,1,1]))

# Create random forcing values
n_segment = 313
n_simhr = 120
forcing_vals = np.random.gamma(1,1,n_segment*n_simhr).reshape(n_simhr,n_segment)
gage_vals = np.random.gamma(1,1,1*n_simhr)
gage_times = np.arange(0.0,n_simhr*3600,3600)  # 3600 for converting hr to sec
# input some NaNs to see how persistence module handles missing values
gage_vals[25:73] = np.nan

# pairs of (flowpath ID, usgs station ID)
full_model.set_value('gages', np.array(['135','01012570','160','01013500','324','01012515']))

# Set DA values
# usgs_df
full_model.set_value('usgs_timeslice_discharge', np.array([1,2,3,4,5,6]))
full_model.set_value('usgs_timeslice_stationId', np.array(['01012570', '01013500', '01012515', '01012570', '01013500', '01012515']))
full_model.set_value('usgs_timeslice_time', np.array(['2015-08-16_00:00:00','2015-08-16_00:03:00','2015-08-16_00:03:00',
                                                '2015-08-16_00:15:00', '2015-08-16_00:18:00','2015-08-16_00:18:00']))
full_model.set_value('usgs_timeslice_discharge_quality', np.array([100,100,100,100,100,100]))

# lastobs_df
full_model.set_value('lastobs_discharge', np.array([1,2,3,4,5,6,7,8,9]))
full_model.set_value('lastobs_stationIdInd', np.array([0,0,0,1,1,1,2,2,2])) 
full_model.set_value('lastobs_timeInd', np.array([0,1,2,0,1,2,0,1,2]))
full_model.set_value('lastobs_stationId', np.array(['01012570','01013500','01012515']))
full_model.set_value('lastobs_time', np.array(['2015-08-15_23:15:00','2015-08-15_23:30:00','2015-08-15_23:45:00',
                                               '2015-08-15_23:15:00','2015-08-15_23:30:00','2015-08-15_23:45:00',
                                               '2015-08-15_23:15:00','2015-08-15_23:30:00','2015-08-15_23:45:00']))
full_model.set_value('lastobs_modelTimeAtOutput', np.array(['2015-08-16_00:00:00']))

#use a standalone DA forcing BMI module for creating DA timeseries and related parameters
DAforcing_model = bmi_DAforcing.bmi_DAforcing()
DAforcing_model.initialize(bmi_cfg_file='/home/dongha.kim/github/t-route/test/BMI/bmi_large_example.yaml')
DAforcing_model.set_value('lake_number', np.array([347987]))
DAforcing_model.set_value('reservoir_type', np.array([4]))
DAforcing_model.update()
import pdb; pdb.set_trace()
full_model.set_value('rfc_discharges', np.array( DAforcing_model._model._rfc_timeseries_discharges))


# Dataframes for saving all output
full_fvd = pd.DataFrame()
full_wbdy = pd.DataFrame()

# Loop through hourly timesteps calling update_until
for hr in range(120):
    # Full network
    # Set dynamic values
    full_model.set_value('land_surface_water_source__volume_flow_rate', forcing_vals[hr,:])
    #full_model.set_value('gage_observations', gage_vals)
    #full_model.set_value('gage_time', gage_times)
    nts = 1
    n = nts*300
    import pdb; pdb.set_trace()
    full_model.update_until(n)
    
    test = lastobs_df_output(
        full_model._model._data_assimilation.lastobs_df,
        full_model._model._time_step,
        nts,
        full_model._model._network.t0 - timedelta(seconds=n),
        full_model._model._network.link_gage_df)

    import pdb; pdb.set_trace()
    '''
    flowveldepth = pd.DataFrame(full_model.get_value('fvd_results').reshape(313,36),index=full_model.get_value('fvd_index'))
    full_fvd = pd.concat([full_fvd, flowveldepth], axis=1)
    '''
import pdb; pdb.set_trace()




