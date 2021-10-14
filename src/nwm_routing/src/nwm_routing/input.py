import pathlib

import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_io as nhd_io
import build_tests  # TODO: Determine whether and how to incorporate this into setup.py

# FIXME
ENV_IS_CL = False
if ENV_IS_CL:
    root = pathlib.Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = pathlib.Path("../../").resolve()


def _input_handler_v03(args):

    custom_input_file = args.custom_input_file
    log_parameters = {}
    supernetwork_parameters = None
    waterbody_parameters = {}
    compute_parameters = {}
    forcing_parameters = {}
    restart_parameters = {}
    output_parameters = {}
    parity_parameters = {}
    data_assimilation_parameters = {}
    diffusive_parameters = {}

    if custom_input_file:
        (
            log_parameters,
            supernetwork_parameters,
            waterbody_parameters,
            compute_parameters,
            forcing_parameters,
            restart_parameters,
            diffusive_parameters,
            output_parameters,
            parity_parameters,
            data_assimilation_parameters,
        ) = nhd_io.read_custom_input_new(custom_input_file)
    else:
        print("CLI input no longer supported")
        raise RuntimeError
    check_inputs(
            log_parameters,
            supernetwork_parameters,
            waterbody_parameters,
            compute_parameters,
            forcing_parameters,
            restart_parameters,
            diffusive_parameters,
            output_parameters,
            parity_parameters,
            data_assimilation_parameters
            )
    return (
        log_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        diffusive_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    )


def _input_handler_v02(args):

    # args = _handle_args()

    custom_input_file = args.custom_input_file
    supernetwork_parameters = None
    waterbody_parameters = {}
    forcing_parameters = {}
    restart_parameters = {}
    output_parameters = {}
    run_parameters = {}
    parity_parameters = {}
    data_assimilation_parameters = {}
    diffusive_parameters = {}

    if custom_input_file:
        (
            supernetwork_parameters,
            waterbody_parameters,
            forcing_parameters,
            restart_parameters,
            output_parameters,
            run_parameters,
            parity_parameters,
            data_assimilation_parameters,
            diffusive_parameters,
            coastal_parameters,
        ) = nhd_io.read_custom_input(custom_input_file)

    else:
        run_parameters["assume_short_ts"] = args.assume_short_ts
        run_parameters["return_courant"] = args.return_courant
        run_parameters["parallel_compute_method"] = args.parallel_compute_method
        run_parameters["subnetwork_target_size"] = args.subnetwork_target_size
        run_parameters["cpu_pool"] = args.cpu_pool
        run_parameters["showtiming"] = args.showtiming
        run_parameters["compute_method"] = args.compute_method

        run_parameters["debuglevel"] = debuglevel = -1 * args.debuglevel
        run_parameters["verbose"] = verbose = args.verbose

        output_parameters["csv_output"] = {}
        output_parameters["csv_output"]["csv_output_folder"] = args.csv_output_folder

        test_folder = pathlib.Path(root, "test")
        geo_input_folder = test_folder.joinpath("input", "geo")

        test_case = args.test_case

        if test_case:

            # call test case assemble function
            (
                supernetwork_parameters,
                run_parameters,
                output_parameters,
                restart_parameters,
                forcing_parameters,
                parity_parameters,
            ) = build_tests.build_test_parameters(
                test_case,
                supernetwork_parameters,
                run_parameters,
                output_parameters,
                restart_parameters,
                forcing_parameters,
                parity_parameters,
            )

        else:
            run_parameters["dt"] = args.dt
            run_parameters["nts"] = args.nts
            run_parameters["qts_subdivisions"] = args.qts_subdivisions
            run_parameters["compute_method"] = args.compute_method

            waterbody_parameters[
                "break_network_at_waterbodies"
            ] = args.break_network_at_waterbodies

            data_assimilation_parameters[
                "data_assimilation_parameters_folder"
            ] = args.data_assimilation_parameters_folder
            data_assimilation_parameters[
                "data_assimilation_filter"
            ] = args.data_assimilation_filter
            data_assimilation_parameters[
                "data_assimilation_csv"
            ] = args.data_assimilation_csv

            restart_parameters[
                "wrf_hydro_channel_restart_file"
            ] = args.wrf_hydro_channel_restart_file
            restart_parameters[
                "wrf_hydro_channel_ID_crosswalk_file"
            ] = args.wrf_hydro_channel_ID_crosswalk_file
            restart_parameters[
                "wrf_hydro_channel_ID_crosswalk_file_field_name"
            ] = args.wrf_hydro_channel_ID_crosswalk_file_field_name
            restart_parameters[
                "wrf_hydro_channel_restart_upstream_flow_field_name"
            ] = args.wrf_hydro_channel_restart_upstream_flow_field_name
            restart_parameters[
                "wrf_hydro_channel_restart_downstream_flow_field_name"
            ] = args.wrf_hydro_channel_restart_downstream_flow_field_name
            restart_parameters[
                "wrf_hydro_channel_restart_depth_flow_field_name"
            ] = args.wrf_hydro_channel_restart_depth_flow_field_name

            forcing_parameters["qlat_const"] = float(args.qlat_const)
            forcing_parameters["qlat_input_folder"] = args.qlat_input_folder
            forcing_parameters["qlat_input_file"] = args.qlat_input_file
            forcing_parameters[
                "qlat_file_pattern_filter"
            ] = args.qlat_file_pattern_filter
            forcing_parameters["qlat_file_index_col"] = args.qlat_file_index_col
            forcing_parameters["qlat_file_value_col"] = args.qlat_file_value_col

            supernetwork = args.supernetwork

        # STEP 0.5: Obtain Supernetwork Parameters for test cases
        if not supernetwork_parameters:
            supernetwork_parameters = nnu.set_supernetwork_parameters(
                supernetwork=supernetwork,
                geo_input_folder=geo_input_folder,
                verbose=False,
                debuglevel=debuglevel,
            )

    return (
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        parity_parameters,
        data_assimilation_parameters,
        diffusive_parameters,
        coastal_parameters,
    )


def check_inputs(
            log_parameters,
            supernetwork_parameters,
            waterbody_parameters,
            compute_parameters,
            forcing_parameters,
            restart_parameters,
            diffusive_parameters,
            output_parameters,
            parity_parameters,
            data_assimilation_parameters
            ):

    #MANDATORY
    try:
        log_parameters['debuglevel']
    except KeyError as err:
        errmsg = 'No debuglevel selected. The default value is 0. Please enter a debuglevel inside log_parameters. Example: debuglevel: 1'
        print(errmsg)
        raise err

    try:
        J = pathlib.Path(supernetwork_parameters['geo_file_path'])
        assert J.is_file() == True
    except AssertionError:
        print("Aborting simulation because forcing file", J, "cannot be not found.")
        raise
    except KeyError as err:
        errmsg = 'No geo_file_path selected. Select a string, file path to directory containing channel geometry data. Example: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc"'
        print(errmsg)
        raise err

    try:
        J = pathlib.Path(supernetwork_parameters['geo_file_path'])
        assert J.is_file() == True
    except AssertionError:
        print("Aborting simulation because forcing file", J, "cannot be not found.")
        raise
    except KeyError as err:
        errmsg = 'No geo_file_path selected. Select a string, file path to directory containing channel geometry data. Example: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc"'
        print(errmsg)
        raise err

    columns = ['key', 'downstream', 'dx', 'n', 'ncc', 's0', 'bw', 'waterbody', 'gages', 'tw', 'twcc', 'alt', 'musk', 'musx', 'cs']   
    col_definitions = ['unique segment identifier',
    'unique identifier of downstream segment',
    'segment length',
    'manning roughness of main channel',
    'mannings roughness of compound channel',
    'channel slope',
    'channel bottom width',
    'waterbody identifier',
    'channel top width',
    'compound channel top width',
    'channel bottom altitude',
    'muskingum K parameter',
    'muskingum X parameter',
    'channel sideslope',
    'gage ID']
    col_examples = ['key: "link"',
            'downstream: "to"',
            'dx: "Length"',
            'n: "n"',  
            'ncc: "nCC"',  
            's0: "So"',  
            'bw: "BtmWdth"',  
            'waterbody: "NHDWaterbodyComID"',
            'gages: "gages"',
            'tw: "TopWdth"',  
            'twcc: "TopWdthCC"',  
            'alt: "alt"',
            'musk: "MusK"',
            'musx: "MusX"',
            'cs: "ChSlp"']
    try:
        for loc, col in enumerate(columns):
            supernetwork_parameters['columns'][col]
    except KeyError as err:
        errmsg = '{} not selected. Select a string, {}. Example: {}'.format(col,col_definitions[loc],col_examples[loc])
        print(errmsg)
        raise err


    if supernetwork_parameters['geo_file_path'][-13:] == "Route_Link.nc":
        try:
            supernetwork_parameters['waterbody_null_code']
        except KeyError as err:
            errmsg = 'No waterbody_null_code selected. Select a null value. The coding in channel gemetry dataset for non waterbody segments under attribute named in `columns: waterbody`. Example: -9999'
            print(errmsg)
            raise err

    #How should we check if a run contains waterbodies? Could write a check then to make sure break_network_at_waterbodies is set to True.
    if waterbody_parameters['break_network_at_waterbodies']:
        columns = ['level_pool_waterbody_parameter_file_path', 'level_pool_waterbody_id', 'level_pool_waterbody_area', 'level_pool_weir_elevation',
        'level_pool_waterbody_max_elevation', 'level_pool_outfall_weir_coefficient', 'level_pool_outfall_weir_length', 'level_pool_overall_dam_length',
        'level_pool_orifice_elevation', 'level_pool_orifice_coefficient', 'level_pool_orifice_area']   
        col_definitions = ['filepath to waterbody parameter file (LAKEPARM.nc)',
        'level_pool_waterbody_id',
        'level_pool_waterbody_area',
        'level_pool_weir_elevation',
        'level_pool_waterbody_max_elevation',
        'level_pool_outfall_weir_coefficient',
        'level_pool_outfall_weir_length',
        'level_pool_overall_dam_length',
        'level_pool_orifice_elevation',
        'level_pool_orifice_coefficient',
        'level_pool_orifice_area']
        col_examples = ['level_pool_waterbody_parameter_file_path: tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/LAKEPARM.nc',
                'level_pool_waterbody_id: lake_id',
                'level_pool_waterbody_area: LkArea',
                'level_pool_weir_elevation: WeirE',
                'level_pool_waterbody_max_elevation: LkMxE',
                'level_pool_outfall_weir_coefficient: WeirC',
                'level_pool_outfall_weir_length: WeirL',
                'level_pool_overall_dam_length: DamL',
                'level_pool_orifice_elevation: OrificeE',
                'level_pool_orifice_coefficient: OrificeC',
                'level_pool_orifice_area: OrificeA']
        try:
            for loc, col in enumerate(columns):
                waterbody_parameters['level_pool'][col]
        except KeyError as err:
            errmsg = '{} not selected. Select a string, the {} column name. Example: {}'.format(col,col_definitions[loc],col_examples[loc])
            print(errmsg)
            raise err

    try:
        compute_parameters['compute_kernel']
    except KeyError as err:
        errmsg = 'No compute_kernel selected. Will default to "V02-structured" - Muskingum Cunge. \
        Please enter a compute_kernel inside compute_parameters. Example: compute_kernel: V02-structured'
        print(errmsg)
        options = ['V02-structured - Muskingum Cunge','diffusive - Diffusive with adaptive timestepping' ,'diffusice_cnt - Diffusive with CNT numerical solution']
        for option in options:
            print("Options: {}".format(option))
        raise err

    if compute_parameters['restart_parameters']['wrf_hydro_channel_restart_file']:
        columns = ['wrf_hydro_channel_ID_crosswalk_file', 'wrf_hydro_channel_ID_crosswalk_file_field_name', 'wrf_hydro_channel_restart_upstream_flow_field_name',
        'wrf_hydro_channel_restart_downstream_flow_field_name',
        'wrf_hydro_channel_restart_depth_flow_field_name', 'wrf_hydro_waterbody_ID_crosswalk_file', 'wrf_hydro_waterbody_ID_crosswalk_file_field_name',
        'wrf_hydro_waterbody_crosswalk_filter_file', 'wrf_hydro_waterbody_crosswalk_filter_file_field_name']   
        col_definitions = [
        'channel geometry file',
        'segment IDs in restart file',
        'upstream flow in restart file',
        'downstream flow in restart file',
        'depth in restart file',
        'lake parameter file',
        'waterbody ID',
        'channel geometry file',
        'waterbody IDs in channel geometry file']
        col_examples = ['wrf_hydro_channel_ID_crosswalk_file: tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc',
                'wrf_hydro_channel_ID_crosswalk_file_field_name: link',
                'wrf_hydro_channel_restart_upstream_flow_field_name: qlink1',
                'wrf_hydro_channel_restart_downstream_flow_field_name: qlink2',
                'wrf_hydro_channel_restart_depth_flow_field_name: hlink',
                'wrf_hydro_waterbody_ID_crosswalk_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/LAKEPARM.nc"',
                'wrf_hydro_waterbody_ID_crosswalk_file_field_name: lake_id',
                'wrf_hydro_waterbody_crosswalk_filter_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc"',
                'wrf_hydro_waterbody_crosswalk_filter_file_field_name: NHDWaterbodyComID']
        try:
            for loc, col in enumerate(columns):
                waterbody_parameters['level_pool'][col]
        except KeyError as err:
            errmsg = '{} not selected. Select a string, the {} column name or filepath. Example: {}'.format(col,col_definitions[loc],col_examples[loc])
            print(errmsg)
            raise err
        import pdb; pdb.set_trace()
    # try:
    #     supernetwork_parameters['synthetic_wb_segments']
    # except KeyError as err:
    #     errmsg = 'No synthetic_wb_segments selected. Select a synthetic values if using NWM 2.1 or 3.0 Route_Link.nc. Synthetic waterbody segment IDs that are used to construct the Great Lakes. \
    #     These segments appear in the NWM 2.1 and 3.0 Route_Link.nc files but are not needed in the routing computation. Example: synthetic_wb_segments:'
    #     print(errmsg)
    #     example_list = [4800002,4800004,4800006,4800007]
    #     for ex in example_list:
    #         print("- {}".format(ex))
    #     raise err


        # check that all forcing files exist
# for f in forcing_filename_list:
#     try:
#         J = pathlib.Path(qlat_input_folder.joinpath(f))     
#         assert J.is_file() == True
#     except AssertionError:
#         print("Aborting simulation because forcing file", J, "cannot be not found.")
#         raise

# try:
#     qlat_input_folder = pathlib.Path(qlat_input_folder)
#     assert qlat_input_folder.is_dir() == True
# except TypeError:
#     print("Aborting simulation because no qlat_input_folder is specified in the forcing_parameters section of the .yaml control file.")
#     raise
# except AssertionError:
#     print("Aborting simulation because the qlat_input_folder:", qlat_input_folder,"does not exist. Please check the the qlat_input_folder variable is correctly entered in the .yaml control file")
#     raise


# try:
#     do_something_that_might_raise_an_exception()
# except ValueError as err:
#     errmsg = 'My custom error message.'
#     raise err(errmsg)