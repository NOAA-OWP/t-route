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
    if compute_parameters['yaml_check']:
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
                compute_parameters['restart_parameters'][col]
        except KeyError as err:
            errmsg = '{} not selected. Select a string, the {} column name or filepath. Example: {}'.format(col,col_definitions[loc],col_examples[loc])
            print(errmsg)
            raise err

    try:
        forcing_parameters['qts_subdivisions']
    except KeyError as err:
        errmsg = 'No qts_subdivisions selected. Please enter a qts_subdivisions inside forcing_parameters. Example: qts_subdivisions: 1. \
        If dt_qlateral = 3600 secs, and dt = 300 secs, then qts_subdivisions = 3600/300 = 12'
        print(errmsg)
        raise err

    try:
        forcing_parameters['dt']
    except KeyError as err:
        errmsg = 'No dt selected. Please enter a dt inside forcing_parameters. Routing simulation time interval. \
        this may be the actual timestep of the numerical solution, but will definitely be the timestep at which \
        flow and depth results are returned from the compute kernel. Example: dt: 300. If dt_qlateral = 3600 secs, and dt \
        = 300 secs, then qts_subdivisions = 3600/300 = 12'
        print(errmsg)
        raise err

    try:
        forcing_parameters['qlat_input_folder']
    except KeyError as err:
        errmsg = 'No qlat_input_folder selected. Select a string, file path to directory containing channel forcing data. Defaults to None and zero-valued lateral inflows are used. Example: qlat_input_folder: "tmp_florence_run_nudging/run_nudging_nwm_channel-only"'
        print(errmsg)
        raise err

    try:
        forcing_parameters['qlat_file_index_col']
    except KeyError as err:
        errmsg = 'No qlat_file_index_col selected. Please select the field name of segment ID in qlateral data.  Example: qlat_file_index_col: feature_id'
        print(errmsg)
        raise err

    try:
        forcing_parameters['qlat_file_value_col']
    except KeyError as err:
        errmsg = 'No qlat_file_value_col selected. Please select the field name of lateral inflow in qlateral data.  Example: qlat_file_value_col: q_lateral'
        print(errmsg)
        raise err

    if forcing_parameters['qlat_forcing_sets']:
        for i in range(0,len(forcing_parameters['qlat_forcing_sets']),1): 
            try:
                forcing_parameters['qlat_forcing_sets'][i]['nts']
            except KeyError as err:
                errmsg = 'No nts selected. Please select the number of timesteps in a loop iteration. nts and max_loop_size will determine the qlat, data assimilation, and parity files used. Example: nts: 48'
                print(errmsg)
                raise err
    else:
        try:
            forcing_parameters['nts']
        except KeyError as err:
            errmsg = 'No nts selected. Please select the number of timesteps in a loop iteration. nts and max_loop_size will determine the qlat, data assimilation, and parity files used. Example: nts: 48'
            print(errmsg)
            raise err
            
        try:
            forcing_parameters['max_loop_size']
        except KeyError as err:
            errmsg = 'No max_loop_size selected. Please select the number of timeslice files in a loop iteration.  Example: max_loop_size: 2 '
            print(errmsg)
            raise err


    try:
        data_assimilation_parameters['data_assimilation_timeslices_folder']
    except KeyError as err:
        errmsg = 'No data_assimilation_timeslices_folder selected. Select a string, file path to directory containing timeSlice files. Example: data_assimilation_timeslices_folder: "nudgingTimeSliceObs_calibration"'
        print(errmsg)
        raise err

    try:
        J = pathlib.Path(data_assimilation_parameters['wrf_hydro_da_channel_ID_crosswalk_file'])
        assert J.is_file() == True
    except AssertionError:
        print("Aborting simulation because assimilation file", J, "cannot be not found.")
        raise
    except KeyError as err:
        errmsg = 'No wrf_hydro_da_channel_ID_crosswalk_file selected. Select a string, file path to channel geometry file. Example: wrf_hydro_da_channel_ID_crosswalk_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc"'
        print(errmsg)
        raise err

    try:
        J = pathlib.Path(data_assimilation_parameters['wrf_hydro_lastobs_file'])
        assert J.is_file() == True
    except AssertionError:
        print("Aborting simulation because assimilation file", J, "cannot be not found.")
        raise
    except KeyError as err:
        errmsg = 'No wrf_hydro_lastobs_file selected. Select a string, filepath to lastobs file. Example: wrf_hydro_lastobs_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/nudgingLastObs.2018-09-01_00:00:00.nc"'
        print(errmsg)
        raise err

    try:
        data_assimilation_parameters['wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time']
    except KeyError as err:
        errmsg = 'No wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time selected. Please select lead time of lastobs relative to simulation start time (secs).  Example: wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time: 0'
        print(errmsg)
        raise err


    try:
        data_assimilation_parameters['lastobs_output_folder']
    except KeyError as err:
        errmsg = 'No lastobs_output_folder selected. Select a string, filepath to lastobs output folder. Example: lastobs_output_folder: "output"'
        print(errmsg)
        raise err


    try:
        data_assimilation_parameters['wrf_lastobs_type']
    except KeyError as err:
        errmsg = 'No wrf_lastobs_type selected. Please select lastobs type of simulation.  Example: wrf_lastobs_type: "obs-based"'
        print(errmsg)
        raise err

    if output_parameters['_csv_output']:
        try:
            output_parameters['_csv_output']['csv_output_folder']
        except KeyError as err:
            errmsg = 'No csv_output_folder selected. Select a string, path to directory where csv output will be written. Example: csv_output_folder: "../../test/output/text"'
            print(errmsg)
            raise err

    if output_parameters['_chrtout_output']:
        try:
            output_parameters['_chrtout_output']['wrf_hydro_channel_output_source_folder']
        except KeyError as err:
            errmsg = 'No wrf_hydro_channel_output_source_folder selected. Select a string, path to directory where un-edited CHRTOUT files are located. \
            These are the same files used as forcings Example: wrf_hydro_channel_output_source_folder: "/glade/work/adamw/forcing/florence_event"'
            print(errmsg)
            raise err

        try:
            output_parameters['_chrtout_output']['wrf_hydro_channel_output_source_folder']
        except KeyError as err:
            errmsg = 'No wrf_hydro_channel_output_source_folder selected. Select a string, path to directory where un-edited CHRTOUT files are located. \
            These are the same files used as forcings Example: wrf_hydro_channel_output_source_folder: "../../forcing/florence_event"'
            print(errmsg)
            raise err

    if output_parameters['_hydro_rst_output']:
        try:
            output_parameters['_hydro_rst_output']['wrf_hydro_channel_restart_source_directory']
        except KeyError as err:
            errmsg = 'No wrf_hydro_channel_restart_source_directory selected. Select a string, path to directory where un-edited HYDRO_RST files are located. Example: wrf_hydro_channel_restart_source_directory: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/"'
            print(errmsg)
            raise err


    if output_parameters['wrf_hydro_parity_check']:
        try:
            output_parameters['wrf_hydro_parity_check']['parity_check_input_folder']
        except KeyError as err:
            errmsg = 'No parity_check_input_folder selected. Select a string, path to directory where WRF-Hydro routed flows are stored in CHRTOUT files. Example: parity_check_input_folder: "tmp_florence_run_nudging/run_nudging_nwm_channel-only"'
            print(errmsg)
            raise err

        try:
            output_parameters['wrf_hydro_parity_check']['parity_check_file_index_col']
        except KeyError as err:
            errmsg = 'No parity_check_file_index_col selected. Please select the name of variable containing segment IDs in CHRTOUT data. Example: parity_check_file_index_col: feature_id'
            print(errmsg)
            raise err

        try:
            output_parameters['wrf_hydro_parity_check']['parity_check_file_value_col']
        except KeyError as err:
            errmsg = 'No parity_check_file_value_col selected. Please select the name of variable containing WRF-Hydro flow in CHRTOUT data. Example: parity_check_file_value_col: streamflow'
            print(errmsg)
            raise err

        try:
            output_parameters['wrf_hydro_parity_check']['parity_check_compare_node']
        except KeyError as err:
            errmsg = 'No parity_check_compare_node selected. Please select the segment ID at which to compare flows. Example: parity_check_compare_node: 8778363'
            print(errmsg)
            raise err
        
        print("YAML check complete")

    