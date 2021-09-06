import json
import pathlib
import pandas as pd
from functools import partial
from datetime import datetime, timedelta
import numpy as np
# TODO: Consider nio and nnw as aliases for these modules...
import troute.nhd_io as nhd_io
import troute.nhd_network as nhd_network
import re

def set_supernetwork_parameters(
    supernetwork="", geo_input_folder=None, verbose=True, debuglevel=0
):
    # TODO: consider managing path concatenation outside this function (and lose the os dependency)

    # The following datasets are extracts from the feature datasets available
    # from https://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/web/data_tools/
    # the CONUS_ge5 and Brazos_LowerColorado_ge5 datasets are included
    # in the github test folder

    supernetwork_options = {
        "Pocono_TEST1",
        "Pocono_TEST2",
        "LowerColorado_Conchos_FULL_RES",
        "Brazos_LowerColorado_ge5",
        "Brazos_LowerColorado_FULL_RES",
        "Brazos_LowerColorado_Named_Streams",
        "CONUS_ge5",
        "Mainstems_CONUS",
        "CONUS_Named_Streams",
        "CONUS_FULL_RES_v20",
        "CapeFear_FULL_RES",
        "Florence_FULL_RES",
        "custom",
    }

    if supernetwork not in supernetwork_options:
        print(
            "Note: please call function with supernetworks set to one of the following:"
        )
        for s in supernetwork_options:
            print(f"'{s}'")
        raise ValueError

    elif supernetwork == "Pocono_TEST1":
        return {
            "geo_file_path": pathlib.Path(
                geo_input_folder, "PoconoSampleData1", "PoconoSampleRouteLink1.shp"
            ).resolve(),
            "columns": {
                "key": "link",
                "downstream": "to",
                "dx": "Length",
                "n": "n",
                "ncc": "nCC",
                "s0": "So",
                "bw": "BtmWdth",
                "waterbody": "NHDWaterbo",
                "tw": "TopWdth",
                "twcc": "TopWdthCC",
                "alt": "alt",
                "musk": "MusK",
                "musx": "MusX",
                "cs": "ChSlp",
            },
            "waterbody_null_code": -9999,
            "title_string": "Pocono Test Example",
            "driver_string": "ESRI Shapefile",
            "terminal_code": 0,
            "layer_string": 0,
            "waterbody_parameter_file_type": "Level_Pool",
            "waterbody_parameter_file_path": pathlib.Path(
                geo_input_folder, "NWM_2.1_Sample_Datasets", "LAKEPARM_CONUS.nc"
            ).resolve(),
            "waterbody_parameter_columns": {
                "waterbody_area": "LkArea",  # area of reservoir
                "weir_elevation": "WeirE",
                "waterbody_max_elevation": "LkMxE",
                "outfall_weir_coefficient": "WeirC",
                "outfall_weir_length": "WeirL",
                "overall_dam_length": "DamL",
                "orifice_elevation": "OrificeE",
                "orifice_coefficient": "OrificeC",
                "orifice_area": "OrificeA",
            },
        }
    elif supernetwork == "Pocono_TEST2":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "Pocono Test 2 Example",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "PoconoRouteLink_TEST2_nwm_mc.txt",
                ).resolve(),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

        # return {
        #'geo_file_path' : pathlib.Path(geo_input_folder
    # , r'PoconoSampleData2'
    # , r'PoconoRouteLink_testsamp1_nwm_mc.shp').resolve()
    # , 'key_col' : 18
    # , 'downstream_col' : 23
    # , 'length_col' : 5
    # , 'manningn_col' : 20
    # , 'manningncc_col' : 21
    # , 'slope_col' : 10
    # , 'bottomwidth_col' : 2
    # , 'topwidth_col' : 11
    # , 'topwidthcc_col' : 12
    # , 'MusK_col' : 6
    # , 'MusX_col' : 7
    # , 'ChSlp_col' : 3
    # , 'terminal_code' : 0
    # , 'title_string' : 'Pocono Test 2 Example'
    # , 'driver_string' : 'ESRI Shapefile'
    # , 'layer_string' : 0
    # }

    elif supernetwork == "LowerColorado_Conchos_FULL_RES":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "NHD 2.0 Conchos Basin of the LowerColorado River",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "LowerColorado_Conchos_FULL_RES.txt",
                ).resolve(),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "Brazos_LowerColorado_ge5":
        return {
            "geo_file_path": pathlib.Path(
                geo_input_folder, "Channels", "NHD_BrazosLowerColorado_Channels.shp"
            ).resolve(),
            "columns": {
                "key": "featureID",
                "downstream": "to",
                "waterbody": "NHDWaterbo",
                "dx": "Length",
                "n": "n",
                "s0": "So",
                "bw": "BtmWdth",
                "musk": "MusK",
                "musx": "MusX",
                "cs": "ChSlp",
            },
            "waterbody_null_code": -9999,
            "title_string": "NHD Subset including Brazos + Lower Colorado\nNHD stream orders 5 and greater",
            "driver_string": "ESRI Shapefile",
            "terminal_code": 0,
            "layer_string": 0,
        }

    elif supernetwork == "Brazos_LowerColorado_FULL_RES":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "NHD 2.0 Brazos and LowerColorado Basins",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "Brazos_LowerColorado_FULL_RES.txt",
                ).resolve(),
                "mask_driver_string": r"csv",
                "mask_layer_string": r"",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "Brazos_LowerColorado_Named_Streams":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "NHD 2.0 GNIS labeled streams in the Brazos and LowerColorado Basins",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "Brazos_LowerColorado_Named_Streams.csv",
                ).resolve(),
                "mask_driver_string": r"csv",
                "mask_layer_string": r"",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "CONUS_ge5":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "NHD CONUS Order 5 and Greater",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder, "Channels", "masks", "CONUS_ge5.txt"
                ).resolve(),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "Mainstems_CONUS":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "CONUS 'Mainstems' (Channels below gages and AHPS prediction points)",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder, r"Channels", r"masks", r"conus_Mainstem_links.txt"
                ).resolve(),
                "mask_driver_string": r"csv",
                "mask_layer_string": r"",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

        # return {
        #     "geo_file_path": pathlib.Path(
        #         geo_input_folder, r"Channels", r"conus_routeLink_subset.nc"
        #     ).resolve(),
        #     "key_col": 0,
        #     "downstream_col": 2,
        #     "length_col": 10,
        #     "manningn_col": 11,
        #     "manningncc_col": 20,
        #     "slope_col": 12,
        #     "bottomwidth_col": 14,
        #     "topwidth_col": 22,
        #     "topwidthcc_col": 21,
        #     "waterbody_col": 15,
        #     "waterbody_null_code": -9999,
        #     "MusK_col": 8,
        #     "MusX_col": 9,
        #     "ChSlp_col": 13,
        #     "terminal_code": 0,
        #     "title_string": 'CONUS "Mainstem"',
        #     "driver_string": "NetCDF",
        #     "layer_string": 0,
        # }

    elif supernetwork == "CONUS_Named_Streams":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "CONUS NWM v2.0 only GNIS labeled streams",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "nwm_reaches_conus_v21_wgnis_name.csv",
                ).resolve(),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "CapeFear_FULL_RES":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "Cape Fear River Basin, NC",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder, "Channels", "masks", "CapeFear_FULL_RES.txt",
                ).resolve(),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "Florence_FULL_RES":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "Hurricane Florence Domain, near Durham NC",  # overwrites other title...
                "mask_file_path": pathlib.Path(
                    geo_input_folder, "Channels", "masks", "Florence_FULL_RES.txt",
                ).resolve(),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "CONUS_FULL_RES_v20":

        ROUTELINK = "RouteLink_NHDPLUS"
        ModelVer = "nwm.v2.0.4"
        ext = "nc"
        sep = "."

        return {
            "geo_file_path": pathlib.Path(
                geo_input_folder, "Channels", sep.join([ROUTELINK, ModelVer, ext])
            ).resolve(),
            "data_link": f"https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/{ModelVer}/parm/domain/{ROUTELINK}{sep}{ext}",
            "columns": {
                "key": "link",
                "downstream": "to",
                "dx": "Length",
                "n": "n",
                "ncc": "nCC",
                "s0": "So",
                "bw": "BtmWdth",
                "tw": "TopWdth",
                "twcc": "TopWdthCC",
                "alt": "alt",
                "waterbody": "NHDWaterbodyComID",
                "musk": "MusK",
                "musx": "MusX",
                "cs": "ChSlp",
            },
            "waterbody_parameter_file_type": "Level_Pool",
            "waterbody_parameter_file_path": pathlib.Path(
                geo_input_folder, "NWM_2.1_Sample_Datasets", "LAKEPARM_CONUS.nc"
            ).resolve(),
            "waterbody_parameter_columns": {
                "waterbody_area": "LkArea",
                "weir_elevation": "WeirE",
                "waterbody_max_elevation": "LkMxE",
                "outfall_weir_coefficient": "WeirC",
                "outfall_weir_length": "WeirL",
                "overall_dam_length": "DamL",
                "orifice_elevation": "OrificeE",
                "oriface_coefficient": "OrificeC",
                "oriface_are": "OrifaceA",
            },
            "waterbody_null_code": -9999,
            "title_string": "CONUS Full Resolution NWM v2.0",
            "driver_string": "NetCDF",
            "terminal_code": 0,
            "layer_string": 0,
        }

    elif supernetwork == "custom":
        custominput = pathlib.Path(geo_input_folder).resolve()
        with open(custominput, "r") as json_file:
            return json.load(json_file)
            # TODO: add error trapping for potentially missing files


def build_connections(supernetwork_parameters):
    cols = supernetwork_parameters["columns"]
    terminal_code = supernetwork_parameters.get("terminal_code", 0)
    synthetic_wb_segments = supernetwork_parameters.get("synthetic_wb_segments", None)
    synthetic_wb_id_offset = supernetwork_parameters.get("synthetic_wb_id_offset", 9.99e11)

    param_df = nhd_io.read(pathlib.Path(supernetwork_parameters["geo_file_path"]))

    pd.set_option('display.max_rows', 500)
    print ("param_df after read")
    print (param_df)
    print ("$$$$$$$$$$$$$$$$$$$$$$$")

    #JDM: Might need to parse out the ints from the fp-id here like in
    #https://github.com/NOAA-OWP/t-route/blob/master/src/external_connections/next_gen_network.py#L82

    param_df = param_df[list(cols.values())]
    param_df = param_df.rename(columns=nhd_network.reverse_dict(cols))
    if synthetic_wb_segments:
        # rename the current key column to key32
        key32_d = {"key":"key32"}
        param_df = param_df.rename(columns=key32_d)
        # create a key index that is int64
        # copy the links into the new column
        param_df["key"] = param_df.key32.astype("int64")
        # update the values of the synthetic reservoir segments
        fix_idx = param_df.key.isin(set(synthetic_wb_segments))
        param_df.loc[fix_idx,"key"] = (param_df[fix_idx].key + synthetic_wb_id_offset).astype("int64")

    param_df = param_df.set_index("key")

    print ("param_df after reindex")
    print (param_df)
    print ("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^6")


    #ngen_nexus_id_to_downstream_comid_mapping_dict = {}
    nexus_to_downstream_flowpath_dict = {}

    if "mask_file_path" in supernetwork_parameters:
        data_mask = nhd_io.read_mask(
            pathlib.Path(supernetwork_parameters["mask_file_path"]),
            layer_string=supernetwork_parameters["mask_layer_string"],
        )

        #print ("data_mask")
        #print (data_mask)
        #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@") 


        param_df = param_df.filter(
            data_mask.iloc[:, supernetwork_parameters["mask_key"]], axis=0
        )

    print ("param_df after read1.5")
    print (param_df)
    print ("$$$$$$$$$$$$$$$$$$$$$$$")
    print ("cols!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print (cols)
    print ("******************")

    #JDM: line below new??
    #param_df = param_df.rename(columns=nhd_network.reverse_dict(cols))
    ####################################################

    print ("param_df after read2")
    print (param_df)
    print ("$$$$$$$$$$$$$$$$$$$$$$$")

    #if "ngen_nexus_id_to_downstream_comid_mapping_json" in supernetwork_parameters:
    #    ngen_nexus_id_to_downstream_comid_mapping_dict = nhd_io.read_ngen_nexus_id_to_downstream_comid_mapping(
    #        pathlib.Path(supernetwork_parameters["ngen_nexus_id_to_downstream_comid_mapping_json"])
    #    )


    if "ngen_nexus_file" in supernetwork_parameters:
        nexus_to_downstream_flowpath_dict = nhd_io.read_nexus_file(
            pathlib.Path(supernetwork_parameters["ngen_nexus_file"])
        )

        #print("ngen_nexus_id_to_downstream_comid_mapping_dict")

        #print("ngen_nexus_id_to_downstream_comid_mapping_dict")
        #print(ngen_nexus_id_to_downstream_comid_mapping_dict)

    #print("param_df")
    #print(param_df)

    #JDM: old param_df line below
    #param_df = param_df.rename(columns=reverse_dict(cols))
    # Rename parameter columns to standard names: from route-link names
    #        key: "link"
    #        downstream: "to"
    #        dx: "Length"
    #        n: "n"  # TODO: rename to `manningn`
    #        ncc: "nCC"  # TODO: rename to `mannningncc`
    #        s0: "So"  # TODO: rename to `bedslope`
    #        bw: "BtmWdth"  # TODO: rename to `bottomwidth`
    #        waterbody: "NHDWaterbodyComID"
    #        gages: "gages"
    #        tw: "TopWdth"  # TODO: rename to `topwidth`
    #        twcc: "TopWdthCC"  # TODO: rename to `topwidthcc`
    #        alt: "alt"
    #        musk: "MusK"
    #        musx: "MusX"
    #        cs: "ChSlp"  # TODO: rename to `sideslope`
    param_df = param_df.sort_index()

    # TODO: Do we need this second, identical call to the one above?
    param_df = param_df.rename(columns=nhd_network.reverse_dict(cols))
    
    print ("param_df after rename")
    print (param_df)
    print ("$$$$$$$$$$$$$$$$$$$$$$$")

    wbodies = {}
    if "waterbody" in cols:
        wbodies = build_waterbodies(
            param_df[["waterbody"]], supernetwork_parameters, "waterbody"
        )
        param_df = param_df.drop("waterbody", axis=1)

    gages = {}
    if "gages" in cols:
        gages = nhd_io.build_filtered_gage_df(param_df[["gages"]])
        param_df = param_df.drop("gages", axis=1)

    # There can be an externally determined terminal code -- that's this first value
    terminal_codes = set()
    terminal_codes.add(terminal_code)
    # ... but there may also be off-domain nodes that are not explicitly identified
    # but which are terminal (i.e., off-domain) as a result of a mask or some other
    # an interior domain truncation that results in a
    # otherwise valid node value being pointed to, but which is masked out or
    # being intentionally separated into another domain.
    terminal_codes = terminal_codes | set(
        param_df[~param_df["downstream"].isin(param_df.index)]["downstream"].values
    )
    connections = nhd_network.extract_connections(
        param_df, "downstream", terminal_codes=terminal_codes
    )
    param_df = param_df.drop("downstream", axis=1)

    param_df = param_df.astype("float32")



    #print ("param_df at end of build_connections")
    #print (param_df)

    #print ("param_df.dtypes")
    #print (param_df.dtypes)


    #print("param_df.index")
    #print(param_df.index)

    #print ("@@@@@@!!!!!!!")

    # datasub = data[['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']]
    #return connections, param_df, wbodies, gages, ngen_nexus_id_to_downstream_comid_mapping_dict
    return connections, param_df, wbodies, gages, nexus_to_downstream_flowpath_dict
    #return connections, param_df, wbodies, gages


def build_waterbodies(
    segment_reservoir_df,
    supernetwork_parameters,
    waterbody_crosswalk_column="waterbody",
):
    """
    segment_reservoir_list
    supernetwork_parameters
    waterbody_crosswalk_column
    """
    wbody_conn = nhd_network.extract_waterbody_connections(
        segment_reservoir_df,
        waterbody_crosswalk_column,
        supernetwork_parameters["waterbody_null_code"],
    )

    # TODO: Add function to read LAKEPARM.nc here
    # TODO: return the lakeparam_df

    return wbody_conn


def organize_independent_networks(connections, wbodies=None):

    rconn = nhd_network.reverse_network(connections)
    independent_networks = nhd_network.reachable_network(rconn)
    reaches_bytw = {}
    for tw, net in independent_networks.items():
        if wbodies:
            path_func = partial(
                nhd_network.split_at_waterbodies_and_junctions, set(wbodies), net
            )
        else:
            path_func = partial(nhd_network.split_at_junction, net)

        reaches_bytw[tw] = nhd_network.dfs_decomposition(net, path_func)

    return independent_networks, reaches_bytw, rconn


def build_channel_initial_state(
    restart_parameters, segment_index=pd.Index([])
):

    channel_restart_file = restart_parameters.get("channel_restart_file", None)

    wrf_hydro_channel_restart_file = restart_parameters.get(
        "wrf_hydro_channel_restart_file", None
    )

    if channel_restart_file:
        q0 = nhd_io.get_channel_restart_from_csv(channel_restart_file)

    elif wrf_hydro_channel_restart_file:

        q0 = nhd_io.get_channel_restart_from_wrf_hydro(
            restart_parameters["wrf_hydro_channel_restart_file"],
            restart_parameters["wrf_hydro_channel_ID_crosswalk_file"],
            restart_parameters["wrf_hydro_channel_ID_crosswalk_file_field_name"],
            restart_parameters["wrf_hydro_channel_restart_upstream_flow_field_name"],
            restart_parameters["wrf_hydro_channel_restart_downstream_flow_field_name"],
            restart_parameters["wrf_hydro_channel_restart_depth_flow_field_name"],
        )
    else:
        # Set cold initial state
        # assume to be zero
        # 0, index=connections.keys(), columns=["qu0", "qd0", "h0",], dtype="float32"
        q0 = pd.DataFrame(
            0, index=segment_index, columns=["qu0", "qd0", "h0"], dtype="float32",
        )
    # TODO: If needed for performance improvement consider filtering mask file on read.
    if not segment_index.empty:
        q0 = q0[q0.index.isin(segment_index)]

    return q0


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

def build_qlateral_array(
    forcing_parameters,
    segment_index=pd.Index([]),
    ts_iterator=None,
    #supernetwork_parameters, #adding this for now, might remove later. Just need to read data_mask
    #ngen_nexus_id_to_downstream_comid_mapping_dict=None,
    nexus_to_downstream_flowpath_dict=None,
    file_run_size=None,
):
    # TODO: set default/optional arguments

    #print ("ngen_nexus_id_to_downstream_comid_mapping_dict2")
    #print (ngen_nexus_id_to_downstream_comid_mapping_dict)

    using_nexus_flows = False

    qts_subdivisions = forcing_parameters.get("qts_subdivisions", 1)
    nts = forcing_parameters.get("nts", 1)
    qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
    qlat_input_file = forcing_parameters.get("qlat_input_file", None)
    nexus_input_folder = forcing_parameters.get("nexus_input_folder", None)
    if qlat_input_folder:
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        if "qlat_files" in forcing_parameters:
            qlat_files = forcing_parameters.get("qlat_files")
            qlat_files = [qlat_input_folder.joinpath(f) for f in qlat_files]
        elif "qlat_file_pattern_filter" in forcing_parameters:
            qlat_file_pattern_filter = forcing_parameters.get(
                "qlat_file_pattern_filter", "*CHRT_OUT*"
            )
            qlat_files = sorted(qlat_input_folder.glob(qlat_file_pattern_filter))

        qlat_file_index_col = forcing_parameters.get(
            "qlat_file_index_col", "feature_id"
        )
        qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", "q_lateral")

        qlat_df = nhd_io.get_ql_from_wrf_hydro_mf(
            qlat_files=qlat_files,
            #ts_iterator=ts_iterator,
            #file_run_size=file_run_size,
            index_col=qlat_file_index_col,
            value_col=qlat_file_value_col,
        )

        qlat_df = qlat_df[qlat_df.index.isin(segment_index)]

    # TODO: These four lines seem extraneous
    #    df_length = len(qlat_df.columns)
    #    for x in range(df_length, 144):
    #        qlat_df[str(x)] = 0
    #        qlat_df = qlat_df.astype("float32")

    elif qlat_input_file:
        qlat_df = nhd_io.get_ql_from_csv(qlat_input_file)


    elif nexus_input_folder:

        using_nexus_flows = True

        #print (nexus_input_folder)
        nexus_input_folder = pathlib.Path(nexus_input_folder)

        if "nexus_file_pattern_filter" in forcing_parameters:
            nexus_file_pattern_filter = forcing_parameters.get(
                "nexus_file_pattern_filter", "nex-*"
            )
            nexus_files = nexus_input_folder.glob(nexus_file_pattern_filter)

            #Declare empty dataframe
            #nexuses_flows_df = pd.DataFrame()

            have_read_in_first_nexus_file = False


            for nexus_file in nexus_files:
                #print (nexus_file)

                split_list = str(nexus_file).split("/")

                #print (split_list)

                nexus_file_name = split_list[-1]

                #print (nexus_file_name)

                nexus_file_name_split = re.split('-|_', nexus_file_name)

                #print (nexus_file_name_split)

                nexus_id = int(nexus_file_name_split[1])

                #print (nexus_id)

                nexus_flows = nhd_io.get_nexus_flows_from_csv(nexus_file)

                #print ("!!!!!!===========nexus_flows-------------")
                #print (nexus_flows)
                
                #comid_df = comid_df.set_index(comid_df.columns[0])

                #if nexus_id in ngen_nexus_id_to_downstream_comid_mapping_dict.keys():
                nexus_flows = nexus_flows.set_index(nexus_flows.columns[0])
                #print ("$$$$$===========nexus_flows-------------")
                #print (nexus_flows)
                

                # Drop original integer index column
                #nexus_flows.drop(nexus_flows.columns[[0]], axis=1, inplace=True)
                #print ("===========nexus_flows-------------")
                #print (nexus_flows)


                nexus_flows = nexus_flows.rename(columns={2: nexus_id})

                #print ("nexus_flows renamed")
                #print (nexus_flows)

                nexus_flows_transposed = nexus_flows.transpose()
                #print ("----------nexus_flows_transposed-------------")
                #print (nexus_flows_transposed)

                # Maybe can change logic for initializing dataframe with append
                if not have_read_in_first_nexus_file:
                    have_read_in_first_nexus_file = True

                    #Need to make the date the index and then do a transformation
                    #to have the date as the header and nex id as the index.
                    #Then append or join each following nexus one.
                    #Then map and reduce to DS segment ids

                    nexuses_flows_df = nexus_flows_transposed

                    # Number of Timesteps plus one
                    # The number of columns in Qlat must be equal to or less than the number of routing timesteps
                    nts = len(nexus_flows) + 1

                    nexus_first_id = nexus_id

                    ##print ("-----------------------------------------")
                    ##print ("nexus_flows.iloc[0]")
                    ##print (nexus_flows.iloc[0])
                    ##print (nexus_flows.iloc[:,0])
                    ##print ("!!!!!!-----------------------------------------")



                    ##print ("nexuses_flows_df")
                    ##print (nexuses_flows_df)

                    #nexuses_flows_df = nexuses_flows_df.set_index(nexus_flows.iloc[:,0])

                    ##print ("nexuses_flows_df after reindex to time ++++++++++++++++++")
                    ##print (nexuses_flows_df)



                else:
                    #TODO: Check on copying and duplication of memory on this??
                    nexuses_flows_df = nexuses_flows_df.append(nexus_flows_transposed)


                    # Number of Timesteps plus one
                    # The number of columns in Qlat must be equal to or less than the number of routing timesteps
                    nts_for_row = len(nexus_flows) + 1

                    if nts_for_row != nts:
                        raise ValueError('Nexus input files number of timesteps discrepancy for nexus-id ' 
                        + str(nexus_first_id) + ' with ' + str(nts) + ' timesteps and nexus-id ' 
                        + str(nexus_id) + ' with ' + str(nts_for_row) + ' timesteps.'
                        )



                ##print ("nexus_flows-------------")
                ##print (nexus_flows)

            print ("@@@@@@@@@@@@nexuses_flows_df")
            print (nexuses_flows_df)

            print ("============================================")

            #Map nexus flows to qlaterals
            #ngen_nexus_id_to_downstream_comid_mapping_dict

            #print (ngen_nexus_id_to_downstream_comid_mapping_dict)

            #nexuses_flows_df

            flowpath_list = []

            #for nexus_key, comid_value in ngen_nexus_id_to_downstream_comid_mapping_dict.items():
            for nexus_key, flowpath_value in nexus_to_downstream_flowpath_dict.items():

                ##print ("nexus_key")
                ##print (nexus_key)

                if flowpath_value not in flowpath_list:
                    flowpath_list.append(flowpath_value)
                else:
                    #print ("%%%%%%%%%%%%")
                    #print (comid_value)
                    delme = 1 

                if flowpath_value not in segment_index:
                    print ("Not in segment_index: " + str(flowpath_value))
                    delme = 1


            #print ("**************")
            #print ("aaaaaaaaaaaaaaaa")

            #print (len(comid_list))


            # Might already be sorted?
            #sorting problem???
            #comid_list = comid_list.sort()

            #print ("comid_list")
            #print (comid_list)


            #comid_df = pd.DataFrame(comid_list)
            #comid_df = pd.DataFrame(comid_list).set_index(comid_list)
           
            #comid_df = comid_df.set_index(comid_df.columns[[0]])
            
            #kinda works???
            #comid_df = comid_df.set_index(comid_df.columns[0])
            
            #comid_df = pd.DataFrame(comid_list, index=[i[0] for i in comid_list])
            #comid_df = pd.DataFrame(comid_list, index=[i for i in comid_list])

            ##print ("comid_df")
            ##print (comid_df)


            already_read_first_nexus_values = False

            #for nexus_key, comid_value in ngen_nexus_id_to_downstream_comid_mapping_dict.items():
            for nexus_key, flowpath_value in nexus_to_downstream_flowpath_dict.items():

                #TODO: simplify below to reduce redundancy in code
                if not already_read_first_nexus_values:
                    already_read_first_nexus_values = True

                    qlat_df_single = pd.DataFrame(nexuses_flows_df.loc[int(nexus_key)])
                    #qlat_df = pd.DataFrame(nexuses_flows_df.loc[int(nexus_key)].transpose())


                    qlat_df_single_transpose = qlat_df_single.transpose()


                    qlat_df_single_transpose = qlat_df_single_transpose.rename(index={int(nexus_key): flowpath_value})

                    #comid_df = comid_df.set_index(comid_df.columns[0])
                    #qlat_df_single_transpose = qlat_df_single_transpose.set_index('1')
                    #qlat_df_single_transpose = qlat_df_single_transpose.set_index(qlat_df_single_transpose.columns[0])

                    #print("qlat_df_single_transpose first") 
                    #print(qlat_df_single_transpose) 

                    qlat_df = qlat_df_single_transpose

                    #print ("qlat_df first")
                    #print (qlat_df)
                    #print ("-------------------------------------------")


                else: 

                    qlat_df_single = pd.DataFrame(nexuses_flows_df.loc[int(nexus_key)])
                    #qlat_df = pd.DataFrame(nexuses_flows_df.loc[int(nexus_key)].transpose())

                    qlat_df_single_transpose = qlat_df_single.transpose()

                    qlat_df_single_transpose = qlat_df_single_transpose.rename(index={int(nexus_key): flowpath_value})

                    #qlat_df_single_transpose = qlat_df_single_transpose.set_index('1')

                    #print("qlat_df_single_transpose") 
                    #print(qlat_df_single_transpose) 

                    #Copying df, memory duplicate????
                    qlat_df = qlat_df.append(qlat_df_single_transpose)


                    #qlat_df = pd.merge(qlat_df, qlat_df_single_transpose
                    
                    #qlat_df = qlat_df.join(qlat_df_single_transpose, how='left')
                    #qlat_df = qlat_df.join(qlat_df_single_transpose, how='outer')
                    
                    #qlat_df = qlat_df.merge(qlat_df_single_transpose, how='outer')

                    pd.set_option('display.max_rows', 500)
                    
                    #print ("qlat_df")
                    #print (qlat_df)
                    #print ("-------------------------------------------")




            #Need to sort qlats

            ############
            #segment_index has full network of segments whereas the downstream segs is a subset of that

  
            full_qlat_df_segment = pd.DataFrame(
                0.0,
                index=segment_index,
                columns=range(nts),
                dtype="float32",
            )

            ##print ("full_qlat_df_segment")
            ##print (full_qlat_df_segment)

            #qlat_df = qlat_df.merge(full_qlat_df_segment, how='right')

            #connection_df['comid'] = connection_df.apply(lambda x: crosswalk_data['cat-' + str(x.name)]['COMID'], axis
            #Need to zero out the values here
            qlat_df_single_transpose_zeros = qlat_df_single_transpose.apply(lambda x: 0.0, axis=0)
            #qlat_df_single_transpose_zeros = qlat_df_single_transpose.apply(np.zeros, axis=1)

            #qlat_df_single_transpose_zeros = qlat_df_single_transpose_zeros.transpose()

            #print ("qlat_df_single_transpose_zeros")
            #print (qlat_df_single_transpose_zeros)

            #maybe transpose in teh to_frame
            qlat_df_single_transpose_zeros_df = qlat_df_single_transpose_zeros.to_frame()


            qlat_df_single_transpose_zeros_df = qlat_df_single_transpose_zeros_df.transpose()


            #print ("qlat_df_single_transpose_zeros_df")
            #print (qlat_df_single_transpose_zeros_df)

            a_segment_index_list = []
            #print ("segment_indexes")
            for a_segment_index in segment_index:
                ##print (a_segment_index)

                if a_segment_index not in a_segment_index_list:
                    a_segment_index_list.append(a_segment_index)

                else:    
                    #print ("repeat segment in mask")
                    #print (a_segment_index)
                    delme = 1

                #if a_segment_index not in comid_list:
                if a_segment_index not in flowpath_list:
                    #add a qlat_df_single_transpose_zeros to qlat_df with the comid          
                    #print ("not in comid_list")
                    #print (a_segment_index)
                    #Copying df, memory duplicate????
                    #qlat_df_single_transpose = qlat_df_single_transpose.rename(index={int(nexus_key): comid_value})
                    qlat_df_single_transpose_zeros_df_renamed = qlat_df_single_transpose_zeros_df.rename(index={0: a_segment_index})
                    
                    #print("qlat_df_single_transpose_zeros_df_renamed")
                    #print(qlat_df_single_transpose_zeros_df_renamed)

                    qlat_df = qlat_df.append(qlat_df_single_transpose_zeros_df_renamed)
                #print ("#############")


            #print ("^^^^^^^^^^^^^^^^^^^")


    else:
        qlat_const = forcing_parameters.get("qlat_const", 0)
        qlat_df = pd.DataFrame(
            qlat_const,
            index=segment_index,
            columns=range(nts // qts_subdivisions),
            dtype="float32",
        )

    pd.set_option('display.max_rows', 500)
    #print ("qlat_df1")
    #print (qlat_df)

    #print ("nts: " + str(nts))


    # TODO: Make a more sophisticated date-based filter
    max_col = 1 + nts // qts_subdivisions

    #print ("max_col: " + str(max_col))

    #print ("len(qlat_df.columns): " + str(len(qlat_df.columns)))

    if len(qlat_df.columns) > max_col:
        qlat_df.drop(qlat_df.columns[max_col:], axis=1, inplace=True)

    #print ("qlat_df1.5")
    #print (qlat_df)

    if not segment_index.empty and not using_nexus_flows:
        qlat_df = qlat_df[qlat_df.index.isin(segment_index)]

    #print ("qlat_df2")
    #print (qlat_df)


    return qlat_df


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

def build_da_sets(data_assimilation_parameters, run_sets, t0):
    
    data_assimilation_timeslices_folder = data_assimilation_parameters.get(
        "data_assimilation_timeslices_folder",
        None
    )
    da_sets = data_assimilation_parameters.get(
        "data_assimilation_sets",
        None
    )
    
    if da_sets:
        pass
        
    elif not data_assimilation_timeslices_folder:
        da_sets = [{} for _ in run_sets]
        
    else:
        data_assimilation_timeslices_folder = pathlib.Path(data_assimilation_timeslices_folder)
        # the number of timeslice files appended to the front- and back-ends
        # of the TimeSlice file interpolation stack
        timeslice_pad = data_assimilation_parameters.get("timeslice_pad",0)

        # timedelta of TimeSlice data - typically 15 minutes
        dt_timeslice = timedelta(minutes = 15)

        da_sets = []
        for (i, set_dict) in enumerate(run_sets):
            da_sets.append({})

            timestamps = pd.date_range(t0 - dt_timeslice * timeslice_pad,
                                       run_sets[i]['final_timestamp']
                                       + dt_timeslice * timeslice_pad,
                                       freq=dt_timeslice)

            filenames = (timestamps.strftime('%Y-%m-%d_%H:%M:%S') 
                        + '.15min.usgsTimeSlice.ncdf').to_list()
            
            # check that all TimeSlice files in the set actually exist
            drop_list = []
            for f in filenames:
                try:
                    J = pathlib.Path(data_assimilation_timeslices_folder.joinpath(f))     
                    assert J.is_file() == True
                except AssertionError:
                    print("Missing TimeSlice file %s", J)
                    drop_list.append(f)
                    
            if drop_list:
                filenames = [x for x in filenames if x not in drop_list]
            
            da_sets[i]['usgs_timeslice_files'] = filenames

            t0 = run_sets[i]['final_timestamp']

    return da_sets
    
def build_data_assimilation(data_assimilation_parameters, run_parameters):
    lastobs_df, da_parameter_dict = build_data_assimilation_lastobs(data_assimilation_parameters)
    usgs_df = build_data_assimilation_usgs_df(data_assimilation_parameters, run_parameters, lastobs_df.index)
    return usgs_df, lastobs_df, da_parameter_dict


def build_data_assimilation_usgs_df(
    data_assimilation_parameters,
    run_parameters,
    lastobs_index=None,
):
    data_assimilation_csv = data_assimilation_parameters.get(
        "data_assimilation_csv", None
    )
    data_assimilation_folder = data_assimilation_parameters.get(
        "data_assimilation_timeslices_folder", None
    )

    usgs_df = pd.DataFrame()
    if not isinstance(lastobs_index, pd.Index):
        lastobs_index = pd.Index()

    if data_assimilation_csv:
        usgs_df = build_data_assimilation_csv(data_assimilation_parameters)
    elif data_assimilation_folder:
        usgs_df = build_data_assimilation_folder(data_assimilation_parameters, run_parameters)

    if not lastobs_index.empty:
        if not usgs_df.empty and not usgs_df.index.equals(lastobs_index):
            print("USGS Dataframe Index Does Not Match Last Observations Dataframe Index")
            usgs_df = usgs_df.loc[lastobs_index]

    return usgs_df


def build_data_assimilation_lastobs(data_assimilation_parameters):
    # TODO: Fix the Logic here according to the following.

    # If there are any observations for data assimilation, there
    # needs to be a complete set in the first time set or else
    # there must be a "LastObs". If there is a complete set in
    # the first time step, the LastObs is optional. If there are
    # no observations for assimilation, there can be a LastObs
    # with an empty usgs dataframe.

    lastobs_df = pd.DataFrame()
    lastobs_file = data_assimilation_parameters.get("wrf_hydro_lastobs_file", None)
    lastobs_start = data_assimilation_parameters.get(
        "wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0
    )
    lastobs_type = data_assimilation_parameters.get("wrf_lastobs_type", "error-based")
    lastobs_crosswalk_file = data_assimilation_parameters.get(
        "wrf_hydro_da_channel_ID_crosswalk_file", None
    )

    if lastobs_file:
        lastobs_df = nhd_io.build_lastobs_df(
            lastobs_file,
            lastobs_crosswalk_file,
            lastobs_type,  # TODO: Confirm that we are using this; delete it if not.
            lastobs_start,
        )

    da_parameter_dict = {}
    da_parameter_dict["da_decay_coefficient"] = data_assimilation_parameters.get("da_decay_coefficient", 120)
    # TODO: Add parameters here for interpolation length (14/59), QC threshold (1.0)

    return lastobs_df, da_parameter_dict


def build_data_assimilation_csv(data_assimilation_parameters):

    usgs_df = nhd_io.get_usgs_df_from_csv(
        data_assimilation_parameters["data_assimilation_csv"],
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
    )

    return usgs_df


def build_data_assimilation_folder(data_assimilation_parameters, run_parameters):

    usgs_timeslices_folder = pathlib.Path(
        data_assimilation_parameters["data_assimilation_timeslices_folder"],
    ).resolve()
    if "data_assimilation_filter" in data_assimilation_parameters:
        da_glob_filter = data_assimilation_parameters["data_assimilation_filter"]
        usgs_files = sorted(usgs_timeslices_folder.glob(da_glob_filter))
    elif "usgs_timeslice_files" in data_assimilation_parameters:
        usgs_files = data_assimilation_parameters.get("usgs_timeslice_files", None)
        usgs_files = [usgs_timeslices_folder.joinpath(f) for f in usgs_files]
    else:
        print("No Files Found for DA")
        # TODO: Handle this with a real exception

    usgs_df = nhd_io.get_usgs_from_time_slices_folder(
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
        usgs_files,
        data_assimilation_parameters.get("qc_threshold", 1),
        data_assimilation_parameters.get("data_assimilation_interpolation_limit", 59),
        run_parameters["dt"],
        run_parameters["t0"],
    )

    return usgs_df
