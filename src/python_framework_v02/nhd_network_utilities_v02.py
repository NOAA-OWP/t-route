import json
import os

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
            "geo_file_path": os.path.join(
                geo_input_folder, "PoconoSampleData1", "PoconoSampleRouteLink1.shp"
            ),
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
            "waterbody_parameter_file_path": os.path.join(
                    geo_input_folder, "NWM_2.1_Sample_Datasets", "LAKEPARM_CONUS.nc"
                ),
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
            }
        }
    elif supernetwork == "Pocono_TEST2":
        rv = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        rv.update(
            {
                "title_string": "Pocono Test 2 Example",  # overwrites other title...
                "mask_file_path": os.path.join(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "PoconoRouteLink_TEST2_nwm_mc.txt",
                ),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

        # return {
        #'geo_file_path' : os.path.join(geo_input_folder
    # , r'PoconoSampleData2'
    # , r'PoconoRouteLink_testsamp1_nwm_mc.shp')
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
                "mask_file_path": os.path.join(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "LowerColorado_Conchos_FULL_RES.txt",
                ),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "Brazos_LowerColorado_ge5":
        return {
            "geo_file_path": os.path.join(
                geo_input_folder, "Channels", "NHD_BrazosLowerColorado_Channels.shp"
            ),
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
                "cs": "ChSlp"
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
                "mask_file_path": os.path.join(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "Brazos_LowerColorado_FULL_RES.txt",
                ),
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
                "mask_file_path": os.path.join(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "Brazos_LowerColorado_Named_Streams.csv",
                ),
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
                "mask_file_path": os.path.join(
                    geo_input_folder, "Channels", "masks", "CONUS_ge5.txt"
                ),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "Mainstems_CONUS":
        dict = set_supernetwork_parameters(
            supernetwork="CONUS_FULL_RES_v20", geo_input_folder=geo_input_folder
        )
        dict.update(
            {
                "title_string": "CONUS 'Mainstems' (Channels below gages and AHPS prediction points)",  # overwrites other title...
                "mask_file_path": os.path.join(
                    geo_input_folder, r"Channels", r"masks", r"conus_Mainstem_links.txt"
                ),
                "mask_driver_string": r"csv",
                "mask_layer_string": r"",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return dict

        # return {
        #     "geo_file_path": os.path.join(
        #         geo_input_folder, r"Channels", r"conus_routeLink_subset.nc"
        #     ),
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
                "mask_file_path": os.path.join(
                    geo_input_folder,
                    "Channels",
                    "masks",
                    "nwm_reaches_conus_v21_wgnis_name.csv",
                ),
                "mask_driver_string": "csv",
                "mask_layer_string": "",
                "mask_key": 0,
                "mask_name": 1,  # TODO: Not used yet.
            }
        )
        return rv

    elif supernetwork == "CONUS_FULL_RES_v20":

        ROUTELINK = "RouteLink_NHDPLUS"
        ModelVer = "nwm.v2.0.2"
        ext = "nc"
        sep = "."

        return {
            "geo_file_path": os.path.join(
                geo_input_folder, "Channels", sep.join([ROUTELINK, ModelVer, ext])
            ),
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
                "waterbody": "NHDWaterbodyComID",
                "musk": "MusK",
                "musx": "MusX",
                "cs": "ChSlp"
            },
            "waterbody_parameter_file_type": "Level_Pool",
            "waterbody_parameter_file_path": os.path.join(
                    geo_input_folder, "NWM_2.1_Sample_Datasets", "LAKEPARM_CONUS.nc"
                ),
            "waterbody_parameter_columns": {
                "waterbody_area": "LkArea",
                "weir_elevation": "WeirE",
                "waterbody_max_elevation": "LkMxE",
                "outfall_weir_coefficient": "WeirC",
                "outfall_weir_length": "WeirL",
                "overall_dam_length": "DamL",
                "orifice_elevation": "OrificeE",
                "oriface_coefficient": "OrificeC",
                "oriface_are": "OrifaceA"
            },
            "waterbody_null_code": -9999,
            "title_string": "CONUS Full Resolution NWM v2.0",
            "driver_string": "NetCDF",
            "terminal_code": 0,
            "layer_string": 0,
        }

    elif supernetwork == "custom":
        custominput = os.path.join(geo_input_folder)
        with open(custominput, "r") as json_file:
            return json.load(json_file)
            # TODO: add error trapping for potentially missing files


def reverse_dict(d):
    """
    Reverse a 1-1 mapping
    Values must be hashable!
    """
    return {v: k for k, v in d.items()}

def get_nhd_connections(supernetwork_parameters={}, verbose=False, debuglevel=0):
    # TODO: convert to get.
    # as in: text=supernetwork_parameters.get("cols_as_text",None):
    # Will need to check if something depends on the None elsewhere.
    if "waterbody_col" not in supernetwork_parameters:
        supernetwork_parameters.update({"waterbody_col": None})
        supernetwork_parameters.update({"waterbody_null_code": None})

    if "mask_file_path" not in supernetwork_parameters:
        # TODO: doing things this way may mean we are reading the same file twice -- fix this [maybe] by implementing an overloaded return
        supernetwork_parameters.update({"mask_file_path": None})
        supernetwork_parameters.update({"mask_layer_string": None})
        supernetwork_parameters.update({"mask_driver_string": None})
        supernetwork_parameters.update({"mask_key_col": None})

    if "data_link" not in supernetwork_parameters:
        supernetwork_parameters.update({"data_link": None})

    return do_connections(
        geo_file_path=supernetwork_parameters["geo_file_path"],
        data_link=supernetwork_parameters["data_link"],
        cols_as_text=supernetwork_parameters.get("cols_as_text", None),
        key_col=supernetwork_parameters["key_col"],
        downstream_col=supernetwork_parameters["downstream_col"],
        length_col=supernetwork_parameters["length_col"],
        terminal_code=supernetwork_parameters["terminal_code"],
        waterbody_col=supernetwork_parameters["waterbody_col"],
        waterbody_null_code=supernetwork_parameters["waterbody_null_code"],
        title_string=supernetwork_parameters["title_string"],
        driver_string=supernetwork_parameters["driver_string"],
        layer_string=supernetwork_parameters["layer_string"],
        mask_file_path=supernetwork_parameters["mask_file_path"],
        mask_layer_string=supernetwork_parameters["mask_layer_string"],
        mask_driver_string=supernetwork_parameters["mask_driver_string"],
        mask_key_col=supernetwork_parameters["mask_key_col"],
        debuglevel=debuglevel,
        verbose=verbose,
    )

def do_connections(
    geo_file_path=None,
    data_link=None,
    cols_as_text=False,
    title_string=None,
    layer_string=None,
    driver_string=None,
    key_col=None,
    downstream_col=None,
    length_col=None,
    terminal_code=None,
    waterbody_col=None,
    waterbody_null_code=None,
    mask_file_path=None,
    mask_driver_string=None,
    mask_layer_string=None,
    mask_key_col=None,
    verbose=False,
    debuglevel=0,
):

    if verbose:
        print(title_string)
    geo_file_rows, geo_keys = get_geo_file_table_rows(
        geo_file_path=geo_file_path,
        data_link=data_link,
        layer_string=layer_string,
        driver_string=driver_string,
        verbose=verbose,
        debuglevel=debuglevel,
    )

    if cols_as_text:
        (
            key_col,
            downstream_col,
            length_col,
            waterbody_col,
        ) = convert_text_cols_to_val(
            geo_keys,
            key_col,
            downstream_col,
            length_col,
            waterbody_col,
            verbose,
            debuglevel,
        )

    if debuglevel <= -1:
        print(f"MASK: {mask_file_path}")
    if mask_file_path:
        mask_file_rows, _ = get_geo_file_table_rows(
            geo_file_path=mask_file_path,
            layer_string=mask_layer_string,
            driver_string=mask_driver_string,
            verbose=verbose,
            debuglevel=debuglevel,
        )
        # TODO: make mask dict with additional attributes, e.g., names
        mask_set = {row[mask_key_col] for row in mask_file_rows}
    else:
        mask_set = {row[key_col] for row in geo_file_rows}

    return build_connections_object(
        geo_file_rows=geo_file_rows,
        mask_set=mask_set,
        key_col=key_col,
        downstream_col=downstream_col,
        length_col=length_col,
        terminal_code=terminal_code,
        waterbody_col=waterbody_col,
        waterbody_null_code=waterbody_null_code,
        verbose=verbose,
        debuglevel=debuglevel,
    )


def set_supernetwork_data(supernetwork="", geo_input_folder=None, verbose=True, debuglevel=0):

    supernetwork_parameters = set_supernetwork_parameters(
        supernetwork=supernetwork, geo_input_folder=geo_input_folder
    )
    supernetwork_values = get_nhd_connections(
        supernetwork_parameters=supernetwork_parameters,
        verbose=verbose,
        debuglevel=debuglevel,
    )
    return supernetwork_parameters, supernetwork_values
