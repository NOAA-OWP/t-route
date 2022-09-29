import json
import pathlib
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

import troute.nhd_io as nhd_io

def build_run_sets(model):
    
    # Create run_sets: sets of forcing files for each loop
    run_sets = build_forcing_sets(model._forcing_parameters, model._network.t0)

    # Create da_sets: sets of TimeSlice files for each loop
    if "data_assimilation_parameters" in model._compute_parameters:
        da_sets = build_da_sets(model._data_assimilation_parameters, run_sets, model._network.t0)
        
    # Create parity_sets: sets of CHRTOUT files against which to compare t-route flows
    if "wrf_hydro_parity_check" in model._output_parameters:
        parity_sets = build_parity_sets(model._parity_parameters, run_sets)
    else:
        parity_sets = []
    
    return run_sets, da_sets, parity_sets

def build_forcings(network, run, forcing_parameters, hybrid_parameters, compute_parameters):
    
    cpu_pool = compute_parameters.get('cpu_pool', None)
    # Create forcing data within network object for first loop iteration
    network.assemble_forcings(run, forcing_parameters, hybrid_parameters, cpu_pool)
    
    return network

def build_data_assimilation(network, data_assimilation_parameters, waterbody_parameters, da_run, forcing_parameters, compute_parameters):
    
    #FIXME: hack to get run_parameters. This is done in input_handler_v2. Probably need
    # to find a better way to do this here though...
    if not 'run_parameters' in locals():
        run_parameters = {'dt': forcing_parameters.get('dt'),
                          'nts': forcing_parameters.get('nts'),
                          'cpu_pool': compute_parameters.get('cpu_pool', None)}
    
    # Create data assimilation object from da_sets for first loop iteration
    data_assimilation = AllDA(data_assimilation_parameters,
                              run_parameters,
                              waterbody_parameters,
                              network,
                              da_run)
    
#    if showtiming:
#        forcing_end_time = time.time()
#        task_times['forcing_time'] += forcing_end_time - network_end_time
    
    return data_assimilation

def build_forcing_sets(
    forcing_parameters,
    t0
):
    
    run_sets = forcing_parameters.get("qlat_forcing_sets", None)
    qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
    nts = forcing_parameters.get("nts", None)
    max_loop_size = forcing_parameters.get("max_loop_size", 12)
    dt = forcing_parameters.get("dt", None)

    try:
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        assert qlat_input_folder.is_dir() == True
    except TypeError:
        raise TypeError("Aborting simulation because no qlat_input_folder is specified in the forcing_parameters section of the .yaml control file.") from None
    except AssertionError:
        raise AssertionError("Aborting simulation because the qlat_input_folder:", qlat_input_folder,"does not exist. Please check the the qlat_input_folder variable is correctly entered in the .yaml control file") from None
    
    forcing_glob_filter = forcing_parameters.get("qlat_file_pattern_filter", "*.CHRTOUT_DOMAIN1")
        
    # TODO: Throw errors if insufficient input data are available

    if run_sets:
        
        # append final_timestamp variable to each set_list
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        for (s, _) in enumerate(run_sets):
            final_chrtout = qlat_input_folder.joinpath(run_sets[s]['qlat_files'
                    ][-1])
            final_timestamp_str = nhd_io.get_param_str(final_chrtout,
                    'model_output_valid_time')
            run_sets[s]['final_timestamp'] = \
                datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')
            
    elif qlat_input_folder:
        
        # Construct run_set dictionary from user-specified parameters
    
        # get the first and seconded files from an ordered list of all forcing files
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        all_files = sorted(qlat_input_folder.glob(forcing_glob_filter))
        first_file = all_files[0]
        second_file = all_files[1]

        # Deduce the timeinterval of the forcing data from the output timestamps of the first
        # two ordered CHRTOUT files
        t1 = nhd_io.get_param_str(first_file, "model_output_valid_time")
        t1 = datetime.strptime(t1, "%Y-%m-%d_%H:%M:%S")
        t2 = nhd_io.get_param_str(second_file, "model_output_valid_time")
        t2 = datetime.strptime(t2, "%Y-%m-%d_%H:%M:%S")
        dt_qlat_timedelta = t2 - t1
        dt_qlat = dt_qlat_timedelta.seconds

        # determine qts_subdivisions
        qts_subdivisions = dt_qlat / dt
        if dt_qlat % dt == 0:
            qts_subdivisions = dt_qlat / dt

        # the number of files required for the simulation
        nfiles = int(np.ceil(nts / qts_subdivisions))

        # list of forcing file datetimes
        datetime_list = [t0 + dt_qlat_timedelta * (n + 1) for n in
                         range(nfiles)]
        datetime_list_str = [datetime.strftime(d, '%Y%m%d%H%M') for d in
                             datetime_list]

        # list of forcing files
        forcing_filename_list = [d_str + ".CHRTOUT_DOMAIN1" for d_str in
                                 datetime_list_str]
        
        # check that all forcing files exist
        for f in forcing_filename_list:
            try:
                J = pathlib.Path(qlat_input_folder.joinpath(f))     
                assert J.is_file() == True
            except AssertionError:
                raise AssertionError("Aborting simulation because forcing file", J, "cannot be not found.") from None
                
        # build run sets list
        run_sets = []
        k = 0
        j = 0
        nts_accum = 0
        nts_last = 0
        while k < len(forcing_filename_list):
            run_sets.append({})

            if k + max_loop_size < len(forcing_filename_list):
                run_sets[j]['qlat_files'] = forcing_filename_list[k:k
                    + max_loop_size]
            else:
                run_sets[j]['qlat_files'] = forcing_filename_list[k:]

            nts_accum += len(run_sets[j]['qlat_files']) * qts_subdivisions
            if nts_accum <= nts:
                run_sets[j]['nts'] = int(len(run_sets[j]['qlat_files'])
                                         * qts_subdivisions)
            else:
                run_sets[j]['nts'] = int(nts - nts_last)

            final_chrtout = qlat_input_folder.joinpath(run_sets[j]['qlat_files'
                    ][-1])
            final_timestamp_str = nhd_io.get_param_str(final_chrtout,
                    'model_output_valid_time')
            run_sets[j]['final_timestamp'] = \
                datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')

            nts_last = nts_accum
            k += max_loop_size
            j += 1
    
    return run_sets

def build_da_sets(da_params, run_sets, t0):
    """
    Create set lists of USGS and/or USACE TimeSlice files used for 
    streamflow and reservoir data assimilation
    
    Arguments
    --------
    - da_params (dict): user-input data assimilation parameters
    - run_sets (list) : forcing files for each run set in the simlation
    - t0 (datetime)   : model initialization time
    
    Returns
    -------
    - da_sets (list)  : lists of USGS and USACE TimeSlice files for each run set 
                        in the simulation
    
    Notes
    -----
    
    """
    
    # check for user-input usgs and usace timeslice directories
    usgs_timeslices_folder = da_params.get(
        "usgs_timeslices_folder",
        None
    )
    usace_timeslices_folder = da_params.get(
        "usace_timeslices_folder",
        None
    )

    # User-specified DA ON/OFF preferences
    usace_da = False
    usgs_da = False
    reservoir_da = da_params.get('reservoir_da', False)
    if reservoir_da:
        usgs_da = reservoir_da.get('reservoir_persistence_usgs', False)
        usace_da = reservoir_da.get('reservoir_persistence_usace', False)
    
    nudging = False
    streamflow_da = da_params.get('streamflow_da', False)
    if streamflow_da:
        nudging = streamflow_da.get('streamflow_nudging', False)
        
    if not usgs_da and not usace_da and not nudging:
        # if all DA capabilities are OFF, return empty dictionary
        da_sets = [{} for _ in run_sets]
    
    # if no user-input timeslice folders, a list of empty dictionaries
    elif not usgs_timeslices_folder and not usace_timeslices_folder:
        # if no timeslice folders, return empty dictionary
        da_sets = [{} for _ in run_sets]
        
    # if user-input timeslice folders are present, build TimeSlice sets
    else:
        
        # create Path objects for each TimeSlice directory
        if usgs_timeslices_folder:
            usgs_timeslices_folder = pathlib.Path(usgs_timeslices_folder)
        if usace_timeslices_folder:
            usace_timeslices_folder = pathlib.Path(usace_timeslices_folder)
        
        # the number of timeslice files appended to the front- and back-ends
        # of the TimeSlice file interpolation stack
        pad_hours = da_params.get("timeslice_lookback_hours",0)
        timeslice_pad = pad_hours * 4 # number of 15 minute TimeSlices in the pad

        # timedelta of TimeSlice data - typically 15 minutes
        dt_timeslice = timedelta(minutes = 15)

        da_sets = [] # initialize list to store TimeSlice set lists

        # Loop through each run set and build lists of available TimeSlice files
        for (i, set_dict) in enumerate(run_sets):
            
            # Append an empty dictionary to the loop, which be used to hold
            # lists of USGS and USACE TimeSlice files.
            da_sets.append({})

            # timestamps of TimeSlice files desired for run set i
            timestamps = pd.date_range(
                t0 - dt_timeslice * timeslice_pad,
                run_sets[i]['final_timestamp'] + dt_timeslice * 4,
                freq=dt_timeslice
            )

            # identify available USGS TimeSlices in run set i
            if (usgs_timeslices_folder and nudging) or (usgs_timeslices_folder and usgs_da):
                filenames_usgs = (timestamps.strftime('%Y-%m-%d_%H:%M:%S') 
                            + '.15min.usgsTimeSlice.ncdf').to_list()
                
                # identify available USGS TimeSlices
                filenames_usgs = _check_timeslice_exists(
                    filenames_usgs, 
                    usgs_timeslices_folder
                )
                
                # Add available TimeSlices to da_sets list
                da_sets[i]['usgs_timeslice_files'] = filenames_usgs
                
            # identify available USACE TimeSlices in run set i
            if usace_timeslices_folder and usace_da:
                filenames_usace = (timestamps.strftime('%Y-%m-%d_%H:%M:%S') 
                            + '.15min.usaceTimeSlice.ncdf').to_list()
                
                # identify available USACE TimeSlices
                filenames_usace = _check_timeslice_exists(
                    filenames_usace, 
                    usace_timeslices_folder
                )
                
                # Add available TimeSlices to da_sets list
                da_sets[i]['usace_timeslice_files'] = filenames_usace

            # reset initialization time for loop set i+1
            t0 = run_sets[i]['final_timestamp']
            
    return da_sets

def build_parity_sets(parity_parameters, run_sets):
    
    parity_sets = parity_parameters.get("parity_check_compare_file_sets", None)
    
    if parity_sets:
        pass
    
    else:
    
        parity_sets = []
        for (i, set_dict) in enumerate(run_sets):
            parity_sets.append({})
            parity_sets[i]['validation_files'] = run_sets[i]['qlat_files']
    
    return parity_sets

def _check_timeslice_exists(filenames, timeslices_folder):
    """
    Check that each TimeSlice file in a list of files exists. Return a list of 
    available files.
    
    Arguments
    ---------
    - filenames               (list of str): TimeSlice filenames
    - timeslices_folder (pathlib.PosixPath): TimeSlice directory
    
    Returns
    -------
    - filenames_existing (list of chr): Existing TimeSlice filenames in 
                                        TimeSlice directory
    
    Notes
    -----
    - This is a utility function used by build_da_sets
    
    """
    
    # check that all USGS TimeSlice files in the set actually exist
    drop_list = []
    for f in filenames:
        try:
            J = pathlib.Path(timeslices_folder.joinpath(f))     
            assert J.is_file() == True
        except AssertionError:
            #LOG.warning("Missing TimeSlice file %s", J)
            drop_list.append(f)

    # Assemble a list of existing TimeSlice files, only
    filenames_existing = [x for x in filenames if x not in drop_list]   
    
    return filenames_existing