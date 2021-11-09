import pathlib
import sys
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_io as nhd_io
import build_tests  # TODO: Determine whether and how to incorporate this into setup.py
from datetime import *
from .log_level_set import log_level_set
import logging


LOG = logging.getLogger('')

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
            preprocessing_parameters,
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
        # TODO: clean _main_.py to remove command line argument comprehension
        LOG.error("CLI input no longer supported")
        raise RuntimeError

    log_level_set(log_parameters)
    LOG = logging.getLogger('')

    if LOG.level <= 10: # DEBUG
        # don't forget to add input checks on user's preprocessing_parameters
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
        preprocessing_parameters,
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
    log_level_set(run_parameters)
    LOG = logging.getLogger('')

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


def _does_file_exist(parameter_name, filepath_input):
    
    debugmsg = 'Checking that a {} parameter is specified and that the file exists'.format(parameter_name)
    LOG.debug(debugmsg)
    try:
        J = pathlib.Path(filepath_input)
        assert J.is_file() == True
    except AssertionError:
        errmsg = '{} parameter entry {} does not exist. Please make sure the file path is correctly entered in the configuration file'.format(parameter_name, filepath_input)
        LOG.error(errmsg)
        raise
    except KeyError as err:
        errmsg = 'A {} parameter must be specified in the configuration file'.format(parameter_name)
        LOG.error(errmsg)
        raise err
    debugmsg  = 'The {} parameter is specified and exists'.format(parameter_name)
    LOG.debug(debugmsg)
    
def _does_path_exist(parameter_name, directory_path_input):
    
    debugmsg = 'Checking that a {} parameter is specified and that the directory exists'.format(parameter_name)
    LOG.debug(debugmsg)
    try:
        J = pathlib.Path(directory_path_input)
        assert J.is_dir() == True
    except AssertionError:
        errmsg = '{} parameter entry {} does not exist. Please make sure the directory path is correctly entered in the configuration file'.format(parameter_name, directory_path_input)
        LOG.error(errmsg)
        raise
    except KeyError as err:
        errmsg = 'A {} parameter must be specified in the configuration file'.format(parameter_name)
        LOG.error(errmsg)
        raise err
    debugmsg  = 'The {} parameter is specified and exists'.format(parameter_name)
    LOG.debug(debugmsg)
    
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

    LOG.debug('***** Begining configuration file (.YAML) check *****')
    
    _does_file_exist('geo_file_path', 
                     supernetwork_parameters['geo_file_path'])
        
    if supernetwork_parameters.get('mask_file_path', None):
        _does_file_exist('mask_file_path', 
                     supernetwork_parameters['mask_file_path'])
    else:
        LOG.debug('No user specified mask file. No mask will be applied to the channel domain.')
        
    if waterbody_parameters.get('break_network_at_waterbodies', None):
        
        # level pool waterbody parameter file
        _does_file_exist('level_pool_waterbody_parameter_file_path', 
                     waterbody_parameters['level_pool']
                         ['level_pool_waterbody_parameter_file_path'])
        
        if waterbody_parameters.get('hybrid_and_rfc', None):
            
            res_usgs = waterbody_parameters['hybrid_and_rfc'].get('reservoir_persistence_usgs', False)
            res_usace = waterbody_parameters['hybrid_and_rfc'].get('reservoir_persistence_usace', False)
            res_forecasts = waterbody_parameters['hybrid_and_rfc'].get('reservoir_rfc_forecasts', False)
            
            if res_usgs or res_usace or res_forecasts:
               
                _does_file_exist('reservoir_parameter_file', 
                         waterbody_parameters['hybrid_and_rfc']
                                 ['reservoir_parameter_file'])
                
                if res_usgs:
                    _does_path_exist('reservoir_usgs_timeslice_path', 
                         waterbody_parameters['hybrid_and_rfc']
                                 ['reservoir_usgs_timeslice_path'])
                    
                if res_usace:
                    _does_path_exist('reservoir_persistence_usace', 
                         waterbody_parameters['hybrid_and_rfc']
                                 ['reservoir_persistence_usace'])
                    
                if res_forecasts:
                    _does_path_exist('reservoir_rfc_forecasts', 
                         waterbody_parameters['hybrid_and_rfc']
                                 ['reservoir_rfc_forecasts'])   
            else:
                LOG.debug('All waterbody DA options are either set to False or not specified. Continuing with no waterbody DA, levelpoool waterbody simulation only.')
                  
    else:
        LOG.debug('break_network_at_waterbodies == False or is not specified, no waterbodies will be simulated')
        
    
    if restart_parameters.get('wrf_hydro_channel_restart_file', None):
        
        _does_file_exist(
            'wrf_hydro_channel_restart_file',
            restart_parameters['wrf_hydro_channel_restart_file']
        ) 
        
        _does_file_exist(
            'wrf_hydro_waterbody_ID_crosswalk_file',
            restart_parameters['wrf_hydro_waterbody_ID_crosswalk_file']
        ) 
        
        _does_file_exist(
            'wrf_hydro_waterbody_crosswalk_filter_file',
            restart_parameters['wrf_hydro_waterbody_crosswalk_filter_file']
        ) 
        
    else:
        LOG.debug('No channel restart file specified, simulation will begin from a cold start')
    
    qlat_folder = forcing_parameters.get('qlat_input_folder',None)
    qlat_sets = forcing_parameters.get('qlat_forcing_sets', None)
    qlat_input_file = forcing_parameters.get('qlat_forcing_sets', None)
    
    if qlat_sets:

        if qlat_sets[0].get('qlat_files', None):
            
            _does_path_exist(
                'qlat_input_folder', 
                forcing_parameters['qlat_input_folder']
            )   
            
            qlat_input_folder = pathlib.Path(forcing_parameters['qlat_input_folder'])
            for (s, _) in enumerate(qlat_sets):
                
                for i, f in enumerate(qlat_sets[s]['qlat_files']):
                    
                    _does_file_exist(
                        'qlat_files[{}][{}]'.format(s, i),
                        qlat_input_folder.joinpath(f)
                    )  
                      
        elif qlat_sets[0].get('qlat_input_file', None):
            
            for (s, _) in enumerate(qlat_sets):
                 
                _does_file_exist(
                    'qlat_input_file',
                    qlat_sets[s]['qlat_input_file']
                )
                
        else:
            errmsg = 'if qlat_forcing_sets variable is in the configuration file, \
            explicitly listed forcing files are expected for each simulation loop. \
            Expected to find qlat_files or qlat_input_file variables in each of the \
            qlat_forcing_sets variable groups.'
            LOG.error(errmsg)
            
    else:
        
        _does_path_exist(
            'qlat_input_folder', 
            forcing_parameters['qlat_input_folder']
        )
    
    if data_assimilation_parameters.get('data_assimilation_timeslices_folder', None):
        
        _does_path_exist(
            'data_assimilation_timeslices_folder', 
            data_assimilation_parameters['data_assimilation_timeslices_folder']
        )
        
        _does_file_exist(
            'wrf_hydro_da_channel_ID_crosswalk_file',
            data_assimilation_parameters['wrf_hydro_da_channel_ID_crosswalk_file']
        )
    
    else:
        LOG.debug('No TimeSlice files proivided for streamflow DA.')
    
    if data_assimilation_parameters.get('wrf_hydro_lastobs_file', None):
        
        _does_file_exist(
            'wrf_hydro_lastobs_file',
            data_assimilation_parameters['wrf_hydro_lastobs_file']
        )
    else:
        LOG.debug('No LastObs file provided for streamflow DA.')
    
    if data_assimilation_parameters.get('lastobs_output_folder', None):
        
        _does_path_exist(
            'lastobs_output_folder', 
            data_assimilation_parameters['lastobs_output_folder']
        )
    else:
        LOG.debug('No LastObs output folder specified. LastObs data will not be written out.')
    
    
    if output_parameters.get('csv_output', None):
        
        if output_parameters['csv_output'].get('csv_output_folder', None):
            
            _does_path_exist(
                'csv_output_folder', 
                output_parameters['csv_output']['csv_output_folder']
            )
            
            output_segs = output_parameters['csv_output'].get('csv_output_segments', None)
            if not output_segs:
                LOG.debug('No csv output segments specified. Results for all domain segments will be written')
            
        else:
            LOG.debug('No csv output folder specified. Results will NOT be written to csv')

    else:
        LOG.debug('No csv output folder specified. Results will NOT be written to csv')
        
    if output_parameters.get('chrtout_output', None):
        
        if output_parameters['chrtout_output'].get('chrtout_read_folder', None):
            
            _does_path_exist(
                'chrtout_read_folder', 
                output_parameters['chrtout_output']['chrtout_read_folder']
            )
            
            if output_parameters['chrtout_output'].get('chrtout_write_folder', None):
                
                _does_path_exist(
                    'chrtout_write_folder', 
                    output_parameters['chrtout_output']['chrtout_write_folder']
                ) 
                
            else:
                LOG.debug('No chrtout_output folder specified. Defaulting CHRTOUT write location to chrtout_read_folder')
                
        else:
            LOG.debug('No chrtout_output parameter specified. Results will NOT be written to CHRTOUT')
    else:
        LOG.debug('No chrtout_output parameter specified. Results will NOT be written to CHRTOUT')
            
    if output_parameters.get('hydro_rst_output', None):
        
        if output_parameters['hydro_rst_output'].get('wrf_hydro_channel_restart_source_directory', None):
            
            _does_path_exist(
                'wrf_hydro_channel_restart_source_directory', 
                output_parameters['hydro_rst_output']['wrf_hydro_channel_restart_source_directory']
            )
            
            if output_parameters['hydro_rst_output'].get('wrf_hydro_channel_restart_output_directory', None):
                
                _does_path_exist(
                    'wrf_hydro_channel_restart_output_directory', 
                    output_parameters['hydro_rst_output']['wrf_hydro_channel_restart_output_directory']
                ) 
                
            else:
                LOG.debug('No wrf_hydro_channel_restart_output_directory folder specified. Defaulting write location to wrf_hydro_channel_restart_source_directory')
                
        else:
            LOG.debug('No hydro_rst_output parameter specified. Restart data will NOT be written out') 
                
    else:
        LOG.debug('No hydro_rst_output parameter specified. Restart data will NOT be written out') 
    
    
    # TODO Parity Check files
    
    # TODO summarize simlation configuration
    
    
    LOG.debug('***** Configuration file check complete *****')


#     columns = ['key', 'downstream', 'dx', 'n', 'ncc', 's0', 'bw', 'waterbody', 'gages', 'tw', 'twcc', 'alt', 'musk', 'musx', 'cs']   
#     col_definitions = ['unique segment identifier',
#     'unique identifier of downstream segment',
#     'segment length',
#     'manning roughness of main channel',
#     'mannings roughness of compound channel',
#     'channel slope',
#     'channel bottom width',
#     'waterbody identifier',
#     'channel top width',
#     'compound channel top width',
#     'channel bottom altitude',
#     'muskingum K parameter',
#     'muskingum X parameter',
#     'channel sideslope',
#     'gage ID']
#     col_examples = ['key: "link"',
#             'downstream: "to"',
#             'dx: "Length"',
#             'n: "n"',  
#             'ncc: "nCC"',  
#             's0: "So"',  
#             'bw: "BtmWdth"',  
#             'waterbody: "NHDWaterbodyComID"',
#             'gages: "gages"',
#             'tw: "TopWdth"',  
#             'twcc: "TopWdthCC"',  
#             'alt: "alt"',
#             'musk: "MusK"',
#             'musx: "MusX"',
#             'cs: "ChSlp"']
    
#     try:
#         for loc, col in enumerate(columns):
#             supernetwork_parameters['columns'][col]
#             # get list of variables from the geofile path
            
#     except KeyError as err:
#         errmsg = '{} ({}) variable name not specified in configuration file. Use the same variable name for {} as it appears in the geo_file_path file. Example: {}'.format(col,col_definitions[loc],col_definitions[loc], col_examples[loc])
#         LOG.error(errmsg)
#         raise err
            
#         columns = ['level_pool_waterbody_parameter_file_path', 'level_pool_waterbody_id', 'level_pool_waterbody_area', 'level_pool_weir_elevation',
#         'level_pool_waterbody_max_elevation', 'level_pool_outfall_weir_coefficient', 'level_pool_outfall_weir_length', 'level_pool_overall_dam_length',
#         'level_pool_orifice_elevation', 'level_pool_orifice_coefficient', 'level_pool_orifice_area']   
#         col_definitions = ['filepath to waterbody parameter file (LAKEPARM.nc)',
#         'level_pool_waterbody_id',
#         'level_pool_waterbody_area',
#         'level_pool_weir_elevation',
#         'level_pool_waterbody_max_elevation',
#         'level_pool_outfall_weir_coefficient',
#         'level_pool_outfall_weir_length',
#         'level_pool_overall_dam_length',
#         'level_pool_orifice_elevation',
#         'level_pool_orifice_coefficient',
#         'level_pool_orifice_area']
#         col_examples = ['level_pool_waterbody_parameter_file_path: tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/LAKEPARM.nc',
#                 'level_pool_waterbody_id: lake_id',
#                 'level_pool_waterbody_area: LkArea',
#                 'level_pool_weir_elevation: WeirE',
#                 'level_pool_waterbody_max_elevation: LkMxE',
#                 'level_pool_outfall_weir_coefficient: WeirC',
#                 'level_pool_outfall_weir_length: WeirL',
#                 'level_pool_overall_dam_length: DamL',
#                 'level_pool_orifice_elevation: OrificeE',
#                 'level_pool_orifice_coefficient: OrificeC',
#                 'level_pool_orifice_area: OrificeA']
#         try:
#             for loc, col in enumerate(columns):
#                 waterbody_parameters['level_pool'][col]
#         except KeyError as err:
#             errmsg = '{} not selected. Select a string, the {} column name. Example: {}'.format(col,col_definitions[loc],col_examples[loc])
#             LOG.error(errmsg)
#             raise err

#     try:
#         compute_parameters['compute_kernel']
#     except KeyError as err:
#         errmsg = 'No compute_kernel selected. Will default to "V02-structured" - Muskingum Cunge. \
#         Please enter a compute_kernel inside compute_parameters. Example: compute_kernel: V02-structured'
#         LOG.error(errmsg)
#         options = ['V02-structured - Muskingum Cunge','diffusive - Diffusive with adaptive timestepping' ,'diffusice_cnt - Diffusive with CNT numerical solution']
#         for option in options:
#             LOG.error("Options: {}".format(option))
#         raise err

#     if compute_parameters['restart_parameters']['wrf_hydro_channel_restart_file']:
#         columns = ['wrf_hydro_channel_ID_crosswalk_file', 'wrf_hydro_channel_ID_crosswalk_file_field_name', 'wrf_hydro_channel_restart_upstream_flow_field_name',
#         'wrf_hydro_channel_restart_downstream_flow_field_name',
#         'wrf_hydro_channel_restart_depth_flow_field_name', 'wrf_hydro_waterbody_ID_crosswalk_file', 'wrf_hydro_waterbody_ID_crosswalk_file_field_name',
#         'wrf_hydro_waterbody_crosswalk_filter_file', 'wrf_hydro_waterbody_crosswalk_filter_file_field_name']   
#         col_definitions = [
#         'channel geometry file',
#         'segment IDs in restart file',
#         'upstream flow in restart file',
#         'downstream flow in restart file',
#         'depth in restart file',
#         'lake parameter file',
#         'waterbody ID',
#         'channel geometry file',
#         'waterbody IDs in channel geometry file']
#         col_examples = ['wrf_hydro_channel_ID_crosswalk_file: tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc',
#                 'wrf_hydro_channel_ID_crosswalk_file_field_name: link',
#                 'wrf_hydro_channel_restart_upstream_flow_field_name: qlink1',
#                 'wrf_hydro_channel_restart_downstream_flow_field_name: qlink2',
#                 'wrf_hydro_channel_restart_depth_flow_field_name: hlink',
#                 'wrf_hydro_waterbody_ID_crosswalk_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/LAKEPARM.nc"',
#                 'wrf_hydro_waterbody_ID_crosswalk_file_field_name: lake_id',
#                 'wrf_hydro_waterbody_crosswalk_filter_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc"',
#                 'wrf_hydro_waterbody_crosswalk_filter_file_field_name: NHDWaterbodyComID']
#         try:
#             for loc, col in enumerate(columns):
#                 compute_parameters['restart_parameters'][col]
#         except KeyError as err:
#             errmsg = '{} not selected. Select a string, the {} column name or filepath. Example: {}'.format(col,col_definitions[loc],col_examples[loc])
#             LOG.error(errmsg)
#             raise err

#     try:
#         forcing_parameters['qts_subdivisions']
#     except KeyError as err:
#         errmsg = 'No qts_subdivisions selected. Please enter a qts_subdivisions inside forcing_parameters. Example: qts_subdivisions: 1. \
#         If dt_qlateral = 3600 secs, and dt = 300 secs, then qts_subdivisions = 3600/300 = 12'
#         LOG.error(errmsg)
#         raise err

#     try:
#         forcing_parameters['dt']
#     except KeyError as err:
#         errmsg = 'No dt selected. Please enter a dt inside forcing_parameters. Routing simulation time interval. \
#         this may be the actual timestep of the numerical solution, but will definitely be the timestep at which \
#         flow and depth results are returned from the compute kernel. Example: dt: 300. If dt_qlateral = 3600 secs, and dt \
#         = 300 secs, then qts_subdivisions = 3600/300 = 12'
#         LOG.error(errmsg)
#         raise err

#     try:
#         forcing_parameters['qlat_input_folder']
#     except KeyError as err:
#         errmsg = 'No qlat_input_folder selected. Select a string, file path to directory containing channel forcing data. Defaults to None and zero-valued lateral inflows are used. Example: qlat_input_folder: "tmp_florence_run_nudging/run_nudging_nwm_channel-only"'
#         LOG.error(errmsg)
#         raise err

#     try:
#         forcing_parameters['qlat_file_index_col']
#     except KeyError as err:
#         errmsg = 'No qlat_file_index_col selected. Please select the field name of segment ID in qlateral data.  Example: qlat_file_index_col: feature_id'
#         LOG.error(errmsg)
#         raise err

#     try:
#         forcing_parameters['qlat_file_value_col']
#     except KeyError as err:
#         errmsg = 'No qlat_file_value_col selected. Please select the field name of lateral inflow in qlateral data.  Example: qlat_file_value_col: q_lateral'
#         LOG.error(errmsg)
#         raise err

#     if forcing_parameters['qlat_forcing_sets']:
#         for i in range(0,len(forcing_parameters['qlat_forcing_sets']),1): 
#             try:
#                 forcing_parameters['qlat_forcing_sets'][i]['nts']
#             except KeyError as err:
#                 errmsg = 'No nts selected. Please select the number of timesteps in a loop iteration. nts and max_loop_size will determine the qlat, data assimilation, and parity files used. Example: nts: 48'
#                 LOG.error(errmsg)
#                 raise err
#     else:
#         try:
#             forcing_parameters['nts']
#         except KeyError as err:
#             errmsg = 'No nts selected. Please select the number of timesteps in a loop iteration. nts and max_loop_size will determine the qlat, data assimilation, and parity files used. Example: nts: 48'
#             LOG.error(errmsg)
#             raise err
            
#         try:
#             forcing_parameters['max_loop_size']
#         except KeyError as err:
#             errmsg = 'No max_loop_size selected. Please select the number of timeslice files in a loop iteration.  Example: max_loop_size: 2 '
#             LOG.error(errmsg)
#             raise err


#     try:
#         data_assimilation_parameters['data_assimilation_timeslices_folder']
#     except KeyError as err:
#         errmsg = 'No data_assimilation_timeslices_folder selected. Select a string, file path to directory containing timeSlice files. Example: data_assimilation_timeslices_folder: "nudgingTimeSliceObs_calibration"'
#         LOG.error(errmsg)
#         raise err

#     try:
#         J = pathlib.Path(data_assimilation_parameters['wrf_hydro_da_channel_ID_crosswalk_file'])
#         assert J.is_file() == True
#     except AssertionError:
#         LOG.error("Aborting simulation because assimilation file", J, "cannot be not found.")
#         raise
#     except KeyError as err:
#         errmsg = 'No wrf_hydro_da_channel_ID_crosswalk_file selected. Select a string, file path to channel geometry file. Example: wrf_hydro_da_channel_ID_crosswalk_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/NWM/DOMAIN/Route_Link.nc"'
#         LOG.error(errmsg)
#         raise err

#     try:
#         J = pathlib.Path(data_assimilation_parameters['wrf_hydro_lastobs_file'])
#         assert J.is_file() == True
#     except AssertionError:
#         LOG.error("Aborting simulation because assimilation file", J, "cannot be not found.")
#         raise
#     except KeyError as err:
#         errmsg = 'No wrf_hydro_lastobs_file selected. Select a string, filepath to lastobs file. Example: wrf_hydro_lastobs_file: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/nudgingLastObs.2018-09-01_00:00:00.nc"'
#         LOG.error(errmsg)
#         raise err

#     try:
#         data_assimilation_parameters['wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time']
#     except KeyError as err:
#         errmsg = 'No wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time selected. Please select lead time of lastobs relative to simulation start time (secs).  Example: wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time: 0'
#         LOG.error(errmsg)
#         raise err


#     try:
#         data_assimilation_parameters['lastobs_output_folder']
#     except KeyError as err:
#         errmsg = 'No lastobs_output_folder selected. Select a string, filepath to lastobs output folder. Example: lastobs_output_folder: "output"'
#         LOG.error(errmsg)
#         raise err


#     try:
#         data_assimilation_parameters['wrf_lastobs_type']
#     except KeyError as err:
#         errmsg = 'No wrf_lastobs_type selected. Please select lastobs type of simulation.  Example: wrf_lastobs_type: "obs-based"'
#         LOG.error(errmsg)
#         raise err

#     if data_assimilation_parameters['wrf_hydro_lastobs_file'] and data_assimilation_parameters['data_assimilation_sets']:
#         data_assimilation_parameters['data_assimilation_sets'][0]['usgs_timeslice_files'][0]
#         lastobs_da_file = datetime.strptime(data_assimilation_parameters['wrf_hydro_lastobs_file'][-22:],'%Y-%m-%d_%H:%M:%S.nc')
#         first_da_file = datetime.strptime(data_assimilation_parameters['data_assimilation_sets'][0]['usgs_timeslice_files'][0][:19],'%Y-%m-%d_%H:%M:%S')
#         if first_da_file != lastobs_da_file:
#             LOG.error("Lastobs file (wrf_hydro_lastobs_file) does not match the first data assimilation file date/time (data_assimilation_sets). Please confirm file dates.")
#         else:
#             pass

#     if output_parameters['csv_output']:
#         try:
#             output_parameters['csv_output']['csv_output_folder']
#         except KeyError as err:
#             errmsg = 'No csv_output_folder selected. Select a string, path to directory where csv output will be written. Example: csv_output_folder: "../../test/output/text"'
#             LOG.error(errmsg)
#             raise err

#     if output_parameters['chrtout_output']:
#         try:
#             output_parameters['chrtout_output']['wrf_hydro_channel_output_source_folder']
#         except KeyError as err:
#             errmsg = 'No wrf_hydro_channel_output_source_folder selected. Select a string, path to directory where un-edited CHRTOUT files are located. \
#             These are the same files used as forcings Example: wrf_hydro_channel_output_source_folder: "/glade/work/adamw/forcing/florence_event"'
#             LOG.error(errmsg)
#             raise err

#         try:
#             output_parameters['chrtout_output']['wrf_hydro_channel_output_source_folder']
#         except KeyError as err:
#             errmsg = 'No wrf_hydro_channel_output_source_folder selected. Select a string, path to directory where un-edited CHRTOUT files are located. \
#             These are the same files used as forcings Example: wrf_hydro_channel_output_source_folder: "../../forcing/florence_event"'
#             LOG.error(errmsg)
#             raise err

#     if output_parameters['hydro_rst_output']:
#         try:
#             output_parameters['hydro_rst_output']['wrf_hydro_channel_restart_source_directory']
#         except KeyError as err:
#             errmsg = 'No wrf_hydro_channel_restart_source_directory selected. Select a string, path to directory where un-edited HYDRO_RST files are located. Example: wrf_hydro_channel_restart_source_directory: "tmp_florence_run_nudging/run_nudging_nwm_channel-only/"'
#             LOG.error(errmsg)
#             raise err


#     if output_parameters['wrf_hydro_parity_check']:
#         try:
#             output_parameters['wrf_hydro_parity_check']['parity_check_input_folder']
#         except KeyError as err:
#             errmsg = 'No parity_check_input_folder selected. Select a string, path to directory where WRF-Hydro routed flows are stored in CHRTOUT files. Example: parity_check_input_folder: "tmp_florence_run_nudging/run_nudging_nwm_channel-only"'
#             LOG.error(errmsg)
#             raise err

#         try:
#             output_parameters['wrf_hydro_parity_check']['parity_check_file_index_col']
#         except KeyError as err:
#             errmsg = 'No parity_check_file_index_col selected. Please select the name of variable containing segment IDs in CHRTOUT data. Example: parity_check_file_index_col: feature_id'
#             LOG.error(errmsg)
#             raise err

#         try:
#             output_parameters['wrf_hydro_parity_check']['parity_check_file_value_col']
#         except KeyError as err:
#             errmsg = 'No parity_check_file_value_col selected. Please select the name of variable containing WRF-Hydro flow in CHRTOUT data. Example: parity_check_file_value_col: streamflow'
#             LOG.error(errmsg)
#             raise err

#         try:
#             output_parameters['wrf_hydro_parity_check']['parity_check_compare_node']
#         except KeyError as err:
#             errmsg = 'No parity_check_compare_node selected. Please select the segment ID at which to compare flows. Example: parity_check_compare_node: 8778363'
#             LOG.error(errmsg)
#             raise err
        
#     LOG.info("YAML check complete")

    