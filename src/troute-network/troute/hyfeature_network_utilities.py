import json
import pathlib
from functools import partial
from datetime import datetime, timedelta
import logging

import pandas as pd
import numpy as np
import netCDF4
from joblib import delayed, Parallel

import troute.nhd_io as nhd_io


LOG = logging.getLogger('')


def build_forcing_sets(
    forcing_parameters,
    t0
):

    run_sets           = forcing_parameters.get("qlat_forcing_sets", None)
    nexus_input_folder = forcing_parameters.get("nexus_input_folder", None)
    nts                = forcing_parameters.get("nts", None)
    max_loop_size      = forcing_parameters.get("max_loop_size", 12)
    dt                 = forcing_parameters.get("dt", None)

    try:
        nexus_input_folder = pathlib.Path(nexus_input_folder)
        assert nexus_input_folder.is_dir() == True
    except TypeError:
        raise TypeError("Aborting simulation because no nexus_input_folder is specified in the forcing_parameters section of the .yaml control file.") from None
    except AssertionError:
        raise AssertionError("Aborting simulation because the nexus_input_folder:", qlat_input_folder,"does not exist. Please check the the nexus_input_folder variable is correctly entered in the .yaml control file") from None

    forcing_glob_filter = forcing_parameters.get("nexus_file_pattern_filter", "*.NEXOUT")
  
    # TODO: Throw errors if insufficient input data are available
    if run_sets:        
        #FIXME: Change it for hyfeature
        '''
        # append final_timestamp variable to each set_list
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        for (s, _) in enumerate(run_sets):
            final_chrtout = qlat_input_folder.joinpath(run_sets[s]['qlat_files'
                    ][-1])
            final_timestamp_str = nhd_io.get_param_str(final_chrtout,
                    'model_output_valid_time')
            run_sets[s]['final_timestamp'] = \
                datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')
        '''  
    elif nexus_input_folder:        
        # Construct run_set dictionary from user-specified parameters

        # get the first and seconded files from an ordered list of all forcing files
        nexus_input_folder = pathlib.Path(nexus_input_folder)
        all_files          = sorted(nexus_input_folder.glob(forcing_glob_filter))
        first_file         = all_files[0]
        second_file        = all_files[1]

        # Deduce the timeinterval of the forcing data from the output timestamps of the first
        # two ordered CHRTOUT files
        df     = pd.read_csv(first_file)
        t1_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
        t1     = datetime.strptime(t1_str,"%Y-%m-%d_%H:%M:%S")
        df     = pd.read_csv(second_file)
        t2_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
        t2     = datetime.strptime(t2_str,"%Y-%m-%d_%H:%M:%S")
        dt_qlat_timedelta = t2 - t1
        dt_qlat = dt_qlat_timedelta.seconds

        # determine qts_subdivisions
        qts_subdivisions = dt_qlat / dt
        if dt_qlat % dt == 0:
            qts_subdivisions = dt_qlat / dt
        # make sure that qts_subdivisions = dt_qlat / dt
        forcing_parameters['qts_subdivisions']= qts_subdivisions

        # the number of files required for the simulation
        nfiles = int(np.ceil(nts / qts_subdivisions))
        
        # list of forcing file datetimes
        #datetime_list = [t0 + dt_qlat_timedelta * (n + 1) for n in
        #                 range(nfiles)]
        # ** Correction ** Because qlat file at time t is constantly applied throughout [t, t+1],
        #               ** n + 1 should be replaced by n
        datetime_list = [t0 + dt_qlat_timedelta * (n) for n in
                         range(nfiles)]        
        datetime_list_str = [datetime.strftime(d, '%Y%m%d%H%M') for d in
                             datetime_list]

        # list of forcing files
        forcing_filename_list = [d_str + "NEXOUT.csv" for d_str in
                                 datetime_list_str]
        
        # check that all forcing files exist
        for f in forcing_filename_list:
            try:
                J = pathlib.Path(nexus_input_folder.joinpath(f))     
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
                run_sets[j]['nexus_files'] = forcing_filename_list[k:k
                    + max_loop_size]
            else:
                run_sets[j]['nexus_files'] = forcing_filename_list[k:]

            nts_accum += len(run_sets[j]['nexus_files']) * qts_subdivisions
            if nts_accum <= nts:
                run_sets[j]['nts'] = int(len(run_sets[j]['nexus_files'])
                                         * qts_subdivisions)
            else:
                run_sets[j]['nts'] = int(nts - nts_last)

            final_nexout        = nexus_input_folder.joinpath(run_sets[j]['nexus_files'
                    ][-1])            
            #final_timestamp_str = nhd_io.get_param_str(final_nexout,
            #        'model_output_valid_time')
            df                  = pd.read_csv(final_nexout)
            final_timestamp_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")           
            
            run_sets[j]['final_timestamp'] = \
                datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')

            nts_last = nts_accum
            k += max_loop_size
            j += 1

    return run_sets

def build_qlateral_array(
    run,
    cpu_pool,
    nexus_to_upstream_flowpath_dict, 
    segment_index=pd.Index([]),
    ts_iterator=None,
    file_run_size=None,
):
    # TODO: set default/optional arguments
    qts_subdivisions = run.get("qts_subdivisions", 1)
    nts = run.get("nts", 1)
    nexus_input_folder = run.get("nexus_input_folder", None)
    qlat_input_file = run.get("qlat_input_file", None)

    if nexus_input_folder:
        nexus_input_folder = pathlib.Path(nexus_input_folder)
        if "nexus_files" in run:
            nexus_files = run.get("nexus_files")
            nexus_files = [nexus_input_folder.joinpath(f) for f in nexus_files]
        elif "nexus_file_pattern_filter" in run:
            nexus_file_pattern_filter = run.get(
                "nexus_file_pattern_filter", "*NEXOUT*"
            )
            nexus_files = sorted(nexus_input_folder.glob(nexus_file_pattern_filter))

        qlat_file_index_col = run.get(
            "qlat_file_index_col", "feature_id"
        )
        qlat_file_value_col = run.get("qlat_file_value_col", "q_lateral")
        gw_bucket_col = run.get("qlat_file_gw_bucket_flux_col","qBucket")
        terrain_ro_col = run.get("qlat_file_terrain_runoff_col","qSfcLatRunoff")


        #nexuses_lateralflows_df = nhd_io.get_ql_from_csv(nexus_files[0])
        '''
        # Parallel reading of qlateral data from CHRTOUT
        with Parallel(n_jobs=cpu_pool) as parallel:
            jobs = []
            for f in qlat_files:
                jobs.append(
                    #delayed(nhd_io.get_ql_from_chrtout)
                    #(f, qlat_file_value_col, gw_bucket_col, terrain_ro_col)
                    delayed(nhd_io.get_ql_from_csv)
                    (f)                    
                )
            ql_list = parallel(jobs)

        # get feature_id from a single CHRTOUT file
        with netCDF4.Dataset(qlat_files[0]) as ds:
            idx = ds.variables[qlat_file_index_col][:].filled()

        # package data into a DataFrame
        qlats_df = pd.DataFrame(
            np.stack(ql_list).T,
            index = idx,
            columns = range(len(qlat_files))
        )    
        '''
        dfs=[]
        for f in nexus_files:
            df = pd.read_csv(f).set_index(['feature_id']) 
            dfs.append(df)
        # lateral flows [m^3/s] are stored at NEXUS points with NEXUS ids
        nexuses_lateralflows_df = pd.concat(dfs, axis=1)  
        
        # Take flowpath ids entering NEXUS and replace NEXUS ids by the upstream flowpath ids 
        qlats_df = pd.concat( (nexuses_lateralflows_df.loc[int(k)].rename(v)
                            for k,v in nexus_to_upstream_flowpath_dict.items() ), axis=1
                            ).T 
        qlats_df.columns=range(len(nexus_files))
        qlats_df = qlats_df[qlats_df.index.isin(segment_index)]
    elif qlat_input_file:
        qlats_df = nhd_io.get_ql_from_csv(qlat_input_file)
    else:
        qlat_const = run.get("qlat_const", 0)
        qlats_df = pd.DataFrame(
            qlat_const,
            index=segment_index,
            columns=range(nts // qts_subdivisions),
            dtype="float32",
        )

    # TODO: Make a more sophisticated date-based filter
    max_col = 1 + nts // qts_subdivisions
    if len(qlats_df.columns) > max_col:
        qlats_df.drop(qlats_df.columns[max_col:], axis=1, inplace=True)

    if not segment_index.empty:
        qlats_df = qlats_df[qlats_df.index.isin(segment_index)]

    return qlats_df