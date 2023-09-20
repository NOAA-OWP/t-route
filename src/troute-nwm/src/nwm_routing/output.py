import time
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
import troute.nhd_io as nhd_io
from build_tests import parity_check
import logging


LOG = logging.getLogger('')

def _reindex_lake_to_link_id(target_df, crosswalk):
    '''
    Utility function for replacing lake ID index values
    with link ID values in a dataframe. This is used to 
    reinedex results dataframes
    
    Arguments:
    ----------
    - target_df (DataFrame): Data frame to be reinexed
    - crosswalk      (dict): Relates lake ids to outlet link ids
    
    Returns:
    --------
    - target_df (DataFrame): Re-indexed with link ids replacing 
                             lake ids
    '''

    # evaluate intersection of lake ids and target_df index values
    # i.e. what are the index positions of lake ids that need replacing?
    lakeids = np.fromiter(crosswalk.keys(), dtype = int)
    idxs = target_df.index.to_numpy()
    lake_index_intersect = np.intersect1d(
        idxs, 
        lakeids, 
        return_indices = True
    )

    # replace lake ids with link IDs in the target_df index array
    linkids = np.fromiter(crosswalk.values(), dtype = int)
    idxs[lake_index_intersect[1]] = linkids[lake_index_intersect[2]]

    # (re) set the target_df index
    target_df.set_index(idxs, inplace = True)
    
    return target_df

def nwm_output_generator(
    run,
    results,
    supernetwork_parameters,
    output_parameters,
    parity_parameters,
    restart_parameters,
    parity_set,
    qts_subdivisions,
    return_courant,
    cpu_pool,
    waterbodies_df,
    waterbody_types_df,
    data_assimilation_parameters=False,
    lastobs_df = None,
    link_gage_df = None,
    link_lake_crosswalk = None,
):

    dt = run.get("dt")
    nts = run.get("nts")
    t0 = run.get("t0")

    if parity_parameters:
        parity_check_waterbody_file = parity_parameters.get("parity_check_waterbody_file", None)
        parity_check_file = parity_parameters.get("parity_check_file", None)
        parity_check_input_folder = parity_parameters.get("parity_check_input_folder", None)
        parity_check_file_index_col = parity_parameters.get(
            "parity_check_file_index_col", None
        )
        parity_check_file_value_col = parity_parameters.get(
            "parity_check_file_value_col", None
        )
        parity_check_compare_node = parity_parameters.get("parity_check_compare_node", None)

        # TODO: find a better way to deal with these defaults and overrides.
        parity_set["parity_check_waterbody_file"] = parity_set.get(
            "parity_check_waterbody_file", parity_check_waterbody_file
        )
        parity_set["parity_check_file"] = parity_set.get(
            "parity_check_file", parity_check_file
        )
        parity_set["parity_check_input_folder"] = parity_set.get(
            "parity_check_input_folder", parity_check_input_folder
        )
        parity_set["parity_check_file_index_col"] = parity_set.get(
            "parity_check_file_index_col", parity_check_file_index_col
        )
        parity_set["parity_check_file_value_col"] = parity_set.get(
            "parity_check_file_value_col", parity_check_file_value_col
        )
        parity_set["parity_check_compare_node"] = parity_set.get(
            "parity_check_compare_node", parity_check_compare_node
        )

    ################### Output Handling
    
    start_time = time.time()

    LOG.info(f"Handling output ...")

    csv_output = output_parameters.get("csv_output", None)
    csv_output_folder = None
    rsrto = output_parameters.get("hydro_rst_output", None)
    chrto = output_parameters.get("chrtout_output", None)
    test = output_parameters.get("test_output", None)
    chano = output_parameters.get("chanobs_output", None)
    wbdyo = output_parameters.get("lakeout_output", None)

    if csv_output:
        csv_output_folder = output_parameters["csv_output"].get(
            "csv_output_folder", None
        )
        csv_output_segments = csv_output.get("csv_output_segments", None)

    
    if csv_output_folder or rsrto or chrto or chano or test or wbdyo:
        
        start = time.time()
        qvd_columns = pd.MultiIndex.from_product(
            [range(nts), ["q", "v", "d"]]
        ).to_flat_index()

        flowveldepth = pd.concat(
            [pd.DataFrame(r[1], index=r[0], columns=qvd_columns) for r in results],
            copy=False,
        )

        if wbdyo and not waterbodies_df.empty:
            
            # create waterbody dataframe for output to netcdf file
            i_columns = pd.MultiIndex.from_product(
                [range(nts), ["i"]]
            ).to_flat_index()

            wbdy = pd.concat(
                [pd.DataFrame(r[6], index=r[0], columns=i_columns) for r in results],
                copy=False,
            )

            wbdy_id_list = waterbodies_df.index.values.tolist()
            flow_df = flowveldepth.loc[wbdy_id_list]
            wbdy = wbdy.loc[wbdy_id_list]

            timestep, variable = zip(*flow_df.columns.tolist())
            timestep_index = np.where(((np.array(list(set(list(timestep)))) + 1) * dt) % (dt * qts_subdivisions) == 0)
            ts_set = set(timestep_index[0].tolist())
            flow_df_col_index = [i for i, e in enumerate(timestep) if e in ts_set]
            flow_df = flow_df.iloc[:,flow_df_col_index]

            timestep, variable = zip(*wbdy.columns.tolist())
            timestep_index = np.where(((np.array(list(set(list(timestep)))) + 1) * dt) % (dt * qts_subdivisions) == 0)
            ts_set = set(timestep_index[0].tolist())
            wbdy_col_index = [i for i, e in enumerate(timestep) if e in ts_set]
            i_df = wbdy.iloc[:,wbdy_col_index]
            q_df = flow_df.iloc[:,0::3]
            d_df = flow_df.iloc[:,2::3]

        # replace waterbody lake_ids with outlet link ids
        if link_lake_crosswalk:
            flowveldepth = _reindex_lake_to_link_id(flowveldepth, link_lake_crosswalk)
            
        # todo: create a unit test by saving FVD array to disk and then checking that
        # it matches FVD array from parent branch or other configurations. 
        # flowveldepth.to_pickle(output_parameters['test_output'])

        if return_courant:
            courant_columns = pd.MultiIndex.from_product(
                [range(nts), ["cn", "ck", "X"]]
            ).to_flat_index()
            courant = pd.concat(
                [
                    pd.DataFrame(r[2], index=r[0], columns=courant_columns)
                    for r in results
                ],
                copy=False,
            )
            
            # replace waterbody lake_ids with outlet link ids
            if link_lake_crosswalk:
                
                # (re) set the flowveldepth index
                courant.set_index(fvdidxs, inplace = True)
            
        LOG.debug("Constructing the FVD DataFrame took %s seconds." % (time.time() - start))

    if test:
        flowveldepth.to_pickle(Path(test))
    
    if wbdyo and not waterbodies_df.empty:
        
        time_index, tmp_variable = map(list,zip(*i_df.columns.tolist()))
        LOG.info("- writing t-route flow results to LAKEOUT files")
        start = time.time()
        for i in range(i_df.shape[1]):            
            nhd_io.write_waterbody_netcdf(
                wbdyo, 
                i_df.iloc[:,[i]],
                q_df.iloc[:,[i]],
                d_df.iloc[:,[i]],
                waterbodies_df,
                waterbody_types_df,
                t0, 
                dt, 
                nts,
                time_index[i],
            )
        
        LOG.debug("writing LAKEOUT files took a total time of %s seconds." % (time.time() - start))
    
    if rsrto:

        LOG.info("- writing restart files")
        start = time.time()
        
        wrf_hydro_restart_dir = rsrto.get(
            "wrf_hydro_channel_restart_source_directory", None
        )

        if wrf_hydro_restart_dir:

            restart_pattern_filter = rsrto.get("wrf_hydro_channel_restart_pattern_filter", "HYDRO_RST.*")
            # list of WRF Hydro restart files
            wrf_hydro_restart_files = sorted(
                Path(wrf_hydro_restart_dir).glob(
                    restart_pattern_filter
                )
            )

            if len(wrf_hydro_restart_files) > 0:
                nhd_io.write_hydro_rst(
                    flowveldepth,
                    wrf_hydro_restart_files,
                    restart_parameters.get("wrf_hydro_channel_restart_file"),
                    dt,
                    nts,
                    t0,
                    restart_parameters.get("wrf_hydro_channel_ID_crosswalk_file"),
                    restart_parameters.get(
                        "wrf_hydro_channel_ID_crosswalk_file_field_name"
                    ),
                )
            else:
                LOG.critical('Did not find any restart files in wrf_hydro_channel_restart_source_directory. Aborting restart write sequence.')
                
        else:
            LOG.critical('wrf_hydro_channel_restart_source_directory not specified in configuration file. Aborting restart write sequence.')
            
        LOG.debug("writing restart files took %s seconds." % (time.time() - start))

    if chrto:
        
        LOG.info("- writing t-route flow results to CHRTOUT files")
        start = time.time()
        
        chrtout_read_folder = chrto.get(
            "wrf_hydro_channel_output_source_folder", None
        )
        
        if chrtout_read_folder:

            chrtout_files = sorted(
                Path(chrtout_read_folder) / f for f in run["qlat_files"]
            )

            nhd_io.write_chrtout(
                flowveldepth,
                chrtout_files,
                qts_subdivisions,
                cpu_pool,
            )
        
        LOG.debug("writing CHRTOUT files took a total time of %s seconds." % (time.time() - start))

    if csv_output_folder: 
    
        LOG.info("- writing flow, velocity, and depth results to .csv")
        start = time.time()

        # create filenames
        # TO DO: create more descriptive filenames
        if supernetwork_parameters.get("title_string", None):
            filename_fvd = (
                "flowveldepth_" + supernetwork_parameters["title_string"] + ".csv"
            )
            filename_courant = (
                "courant_" + supernetwork_parameters["title_string"] + ".csv"
            )
        else:
            run_time_stamp = datetime.now().isoformat()
            filename_fvd = "flowveldepth_" + run_time_stamp + ".csv"
            filename_courant = "courant_" + run_time_stamp + ".csv"

        output_path = Path(csv_output_folder).resolve()

        # no csv_output_segments are specified, then write results for all segments
        if not csv_output_segments:
            csv_output_segments = flowveldepth.index
        
        flowveldepth = flowveldepth.sort_index()
        flowveldepth.loc[csv_output_segments].to_csv(output_path.joinpath(filename_fvd))

        if return_courant:
            courant = courant.sort_index()
            courant.loc[csv_output_segments].to_csv(output_path.joinpath(filename_courant))

        LOG.debug("writing CSV file took %s seconds." % (time.time() - start))
        # usgs_df_filtered = usgs_df[usgs_df.index.isin(csv_output_segments)]
        # usgs_df_filtered.to_csv(output_path.joinpath("usgs_df.csv"))
        
    if chano:
        
        LOG.info("- writing t-route flow results at gage locations to CHANOBS file")
        start = time.time()
        
        # replace waterbody lake_ids with outlet link ids
        if link_lake_crosswalk:
            link_gage_df = _reindex_lake_to_link_id(link_gage_df, link_lake_crosswalk)

        if isinstance(chano['chanobs_output_directory'], Path):
            chano['chanobs_output_directory'] = str(chano['chanobs_output_directory']) + '/'
            chano['chanobs_filepath'] = str(chano['chanobs_filepath'])

        nhd_io.write_chanobs(
            Path(chano['chanobs_output_directory'] + chano['chanobs_filepath']), 
            flowveldepth, 
            link_gage_df, 
            t0, 
            dt, 
            nts,
            # TODO allow user to pass a list of segments at which they would like to print results
            # rather than just printing at gages. 
        )
        
        LOG.debug("writing flow data to CHANOBS took %s seconds." % (time.time() - start))       

        # Write out LastObs as netcdf.
        # This is only needed if 1) streamflow nudging is ON and 2) a lastobs output
        # folder is provided by the user.
        lastobs_output_folder = None
        nudging_true = None
        streamflow_da = data_assimilation_parameters.get('streamflow_da', None)
        if streamflow_da:
            lastobs_output_folder = streamflow_da.get(
                "lastobs_output_folder", None
                )
            nudging_true = streamflow_da.get('streamflow_nudging', None)

        if nudging_true and lastobs_output_folder:

            LOG.info("- writing lastobs files")
            start = time.time()

            nhd_io.lastobs_df_output(
                lastobs_df,
                dt,
                nts,
                t0,
                link_gage_df['gages'],
                lastobs_output_folder,
            )

            LOG.debug("writing lastobs files took %s seconds." % (time.time() - start))

    if 'flowveldepth' in locals():
        LOG.debug(flowveldepth)

    LOG.debug("output complete in %s seconds." % (time.time() - start_time))

    ################### Parity Check

    if parity_set:
        
        LOG.info(
            "conducting parity check, comparing WRF Hydro results against t-route results"
        )
    
        start_time = time.time()
       
        parity_check(
            parity_set, results,
        )

        LOG.debug("parity check complete in %s seconds." % (time.time() - start_time))