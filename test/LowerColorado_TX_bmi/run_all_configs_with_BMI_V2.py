import sys
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle
from datetime import datetime, timedelta
import time
# /home/amin/Amin_fork/t-route/test/Apr2023_testcase/Extended_AnA.yaml
sys.path.append("/home/amin/Amin_fork/t-route/src/")
import bmi_troute
import bmi_DAforcing
import model_DAforcing
# Define generic troute 'run' function for any configuration
def run_troute(config):    
    print('***************',config,' Run***************', sep='')

    base_path = '/home/amin/Amin_fork/t-route/test/LowerColorado_TX_bmi/'
    config_dict = {
        'Extended AnA': {
            'yaml': 'configs/Extended_AnA.yaml',
            'duration': 25,
            'output': 'output/Extended_AnA/bmi/',
            'restart': 23
            },
        'Standard AnA': {
            'yaml': 'configs/Standard_AnA.yaml',
            'duration': 3,
            'output': 'output/Standard_AnA/bmi/',
            'restart': 2
            },
        'ShortRange Forecast': {
            'yaml': 'configs/ShortRange_Forecast.yaml',
            'duration': 18,
            'output': 'output/ShortRange_Forecast/bmi/',
            'restart': None
            },
        'MediumRange Forecast': {
            'yaml': 'configs/MediumRange_Forecast.yaml',
            'duration': 240,
            'output': 'output/MediumRange_Forecast/bmi/',
            'restart': None
            },
        'LongRange Forecast': {
            'yaml': 'configs/LongRange_Forecast.yaml',
            'duration': 720,
            'output': 'output/LongRange_Forecast/bmi/',
            'restart': None
            }
        }
    
    #-------------------------------------
    # Begin our t-route related operations
    #-------------------------------------
    # Initialize DA model with configuration file
    print('Initialize DA Module...')
    DAforcing = bmi_DAforcing.bmi_DAforcing()
    DAforcing.initialize(bmi_cfg_file = base_path + config_dict[config]['yaml'])
    # DAforcing_flat = model_DAforcing.DAforcing_model()
    # Initialize troute model with configuration file
    print('Initialize t-route Module...')
    troute = bmi_troute.bmi_troute()
    troute.initialize(bmi_cfg_file = base_path + config_dict[config]['yaml'])
    
    # Set DA values
    troute.set_value('usgs_df', DAforcing.get_value('usgs_df'))
    troute.set_value('reservoir_usgs_df', DAforcing.get_value('reservoir_usgs_df'))
    troute.set_value('reservoir_usace_df', DAforcing.get_value('reservoir_usace_df'))
    troute.set_value('rfc_timeseries_df', DAforcing.get_value('rfc_timeseries_df'))
    troute.set_value('waterbody_df', DAforcing.get_value('waterbody_df'))
    troute.set_value('waterbody_df_index', DAforcing.get_value('waterbody_df_ids'))
    troute.set_value('lastobs_df', DAforcing.get_value('lastobs_df'))
    troute.set_value('lastobs_df_index', DAforcing.get_value('lastobs_df_ids'))

    # Set restart values
    troute.set_value('q0', DAforcing.get_value('q0'))
    troute.set_value('q0_index', DAforcing.get_value('q0_ids'))
    troute.set_value('t0', DAforcing.get_value('t0'))
    
    # Read forcing data first. This info will eventually be provided by model engine.
    timeseries_start = DAforcing._model._t0
    timeseries_end = timeseries_start + timedelta(hours=config_dict[config]['duration'] - 1)
    delta = timedelta(hours=1)
    forcing_files = []
    while timeseries_start <= timeseries_end:
        forcing_files.append('/home/amin/Amin_fork/t-route/test/LowerColorado_TX_bmi/channel_forcing/' + timeseries_start.strftime('%Y%m%d%H%M') + '.CHRTOUT_DOMAIN1.csv')
        timeseries_start += delta
    
    # Loop through hourly timesteps calling update_until
    print('Begin routing...')
    for hr in range(config_dict[config]['duration']):
        # Retrieve forcing data
        forcing_vals = pd.read_csv(forcing_files[hr]).set_index('feature_id')
        # Get rid of duplicate nexus points:
        idx = forcing_vals.index.duplicated()
        forcing_vals = forcing_vals.loc[~idx]
        
        # Full network
        # Set dynamic values
        troute.set_value('land_surface_water_source__volume_flow_rate', np.array(forcing_vals.iloc[:,0]))
        troute.set_value('land_surface_water_source__id', np.array(forcing_vals.index))

        # Set duration to run model for
        nts = 12
        n = nts*300
        
        # Run our model
        troute.update_until(n)
        if hr!=config_dict[config]['restart']:
            DAforcing.set_value('nudging', troute.get_value('nudging'))
            DAforcing.set_value('nudging_ids', troute.get_value('nudging_ids'))
            DAforcing.set_value('fvd_results', troute.get_value('fvd_results'))
            DAforcing.set_value('fvd_index', troute.get_value('fvd_index'))
            DAforcing.set_value('t-route_model_time', troute.get_current_time())

            DAforcing.update()

           
        if hr==config_dict[config]['restart']:
            DAforcing.set_value('q0', troute.get_value('q0'))
            DAforcing.set_value('q0_ids', troute.get_value('q0_index'))
            
            DAforcing.set_value('waterbody_df', troute.get_value('waterbody_df'))
            DAforcing.set_value('waterbody_df_ids', troute.get_value('waterbody_df_index'))

            DAforcing.set_value('lastobs_df', troute.get_value('lastobs_df'))
            DAforcing.set_value('lastobs_df_ids', troute.get_value('lastobs_df_index'))

            DAforcing.set_value('t-route_model_time', troute.get_current_time())
            DAforcing.set_value('nudging', troute.get_value('nudging'))
            DAforcing.set_value('nudging_ids', troute.get_value('nudging_ids'))
            DAforcing.set_value('fvd_results', troute.get_value('fvd_results'))
            DAforcing.set_value('fvd_index', troute.get_value('fvd_index'))
            DAforcing.set_value('t-route_model_time', troute.get_current_time())
            
            DAforcing.update()

        #print('Extended AnA model time:', troute.get_current_time())
        
    troute._model.print_timing_summary()


##### Run t-route #####
# Time logging
total_start_time = time.time()

# Run 5 configs and save results as pickle files for analysis:
run_troute('Extended AnA')
# run_troute('Standard AnA')
# run_troute('ShortRange Forecast')
# run_troute('MediumRange Forecast')
# run_troute('LongRange Forecast')

total_time = time.time() - total_start_time

print('Total run:', total_time)
print('Total time: {} secs'.format(round(total_time,2)))

