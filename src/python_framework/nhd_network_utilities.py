import networkbuilder as networkbuilder
import os
import traceback
import geopandas as gpd
import pandas as pd
import zipfile
import xarray as xr


def get_geo_file_table_rows(
        geo_file_path = None
        , layer_string = None
        , driver_string = None
        , verbose = False
        , debuglevel = 0
        ):

    # NOTE: these methods can lose the "connections" and "rows" arguments when
    # implemented as class methods where those arrays are members of the class.
    if driver_string == 'NetCDF': # Use Xarray to read a netcdf table
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; driver: {driver_string}')
        try:
            geo_file = xr.open_dataset(geo_file_path)
            geo_file_rows = (geo_file.to_dataframe()).values
        except Exception as e:
            print (e)
            if debuglevel <= -1: traceback.print_exc()
        # The xarray method for NetCDFs was implemented after the geopandas method for 
        # GIS source files. It's possible (probable?) that we are doing something 
        # inefficient by converting away from the Pandas dataframe.
        # TODO: Check the optimal use of the Pandas dataframe
    elif driver_string == 'csv': # Use Pandas to read zipped csv/txt
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; driver: {driver_string}')
        try:
            geo_file = pd.read_csv(geo_file_path)
        except Exception as e:
            print (e)
            if debuglevel <= -1: traceback.print_exc()
        geo_file_rows = geo_file.to_numpy()
    elif driver_string == 'zip': # Use Pandas to read zipped csv/txt
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; driver: {driver_string}')
        try:
            with zipfile.ZipFile(geo_file_path, 'r') as zcsv:
                with zcsv.open(layer_string) as csv:
                    geo_file = pd.read_csv(csv)
        except Exception as e:
            print (e)
            if debuglevel <= -1: traceback.print_exc()
        geo_file_rows = geo_file.to_numpy()
    else: # Read Shapefiles, Geodatabases with Geopandas/fiona
        if debuglevel <= -1: print(f'reading -- dataset: {geo_file_path}; layer: {layer_string}; fiona driver: {driver_string}')
        try:
            geo_file = gpd.read_file(geo_file_path, driver=driver_string, layer=layer_string)
            geo_file_rows = geo_file.to_numpy()
        except Exception as e:
            print (e)
            if debuglevel <= -1: traceback.print_exc()
        if debuglevel <= -2: 
            try: 
                geo_file.plot() 
            except:
                print(r'cannot plot geofile (not necessarily a problem)')
    if debuglevel <= -1: print(geo_file.head()) # Preview the first 5 lines of the loaded data

    return geo_file_rows

def build_connections_object(
        geo_file_rows = None
        , mask_set = None
        , key_col = None
        , downstream_col = None
        , length_col = None
        , terminal_code = None
        , verbose = False
        , debuglevel = 0
        ):
    (connections) = networkbuilder.get_down_connections(
                    rows = geo_file_rows
                    , mask_set = mask_set
                    , key_col = key_col
                    , downstream_col = downstream_col
                    , length_col = length_col
                    , verbose = verbose
                    , debuglevel = debuglevel)
    
    (all_keys, ref_keys, headwater_keys
        , terminal_keys
        , terminal_ref_keys
        , circular_keys) = networkbuilder.determine_keys(
                    connections = connections
                    , key_col = key_col
                    , downstream_col = downstream_col
                    , terminal_code = terminal_code
                    , verbose = verbose
                    , debuglevel = debuglevel)
    
    (junction_keys, visited_keys
     , visited_terminal_keys
     , junction_count) = networkbuilder.get_up_connections(
                    connections = connections
                    , terminal_code = terminal_code
                    , headwater_keys = headwater_keys
                    , terminal_keys = terminal_keys
                    , verbose = verbose
                    , debuglevel = debuglevel)
    return connections, all_keys, ref_keys, headwater_keys \
        , terminal_keys, terminal_ref_keys \
        , circular_keys, junction_keys \
        , visited_keys, visited_terminal_keys \
        , junction_count

def do_connections(
        geo_file_path = None
        , title_string = None
        , layer_string = None
        , driver_string = None
        , key_col = None
        , downstream_col = None
        , length_col = None
        , terminal_code = None
        , mask_file_path = None
        , mask_driver_string = None
        , mask_layer_string = None
        , mask_key_col = None
        , verbose = False
        , debuglevel = 0
        ):

    if verbose: print(title_string)
    geo_file_rows = get_geo_file_table_rows(
        geo_file_path = geo_file_path
        , layer_string = layer_string
        , driver_string = driver_string
        , verbose = verbose
        , debuglevel = debuglevel
    )

    if debuglevel <= -1: print(f'MASK: {mask_file_path}')
    if mask_file_path:
        mask_file_rows = get_geo_file_table_rows(
            geo_file_path = mask_file_path
            , layer_string = mask_layer_string
            , driver_string = mask_driver_string
            , verbose = verbose
            , debuglevel = debuglevel
            )
        #TODO: make mask dict with additional attributes, e.g., names
        mask_set = {row[mask_key_col] for row in mask_file_rows}
    else: mask_set = {row[key_col] for row in geo_file_rows}

    return build_connections_object(
        geo_file_rows = geo_file_rows
        , mask_set = mask_set
        , key_col = key_col
        , downstream_col = downstream_col
        , length_col = length_col
        , terminal_code = terminal_code
        , verbose = verbose
        , debuglevel = debuglevel
        )

    # return connections, all_keys, ref_keys, headwater_keys \
    #     , terminal_keys, terminal_ref_keys \
    #     , circular_keys, junction_keys \
    #     , visited_keys, visited_terminal_keys \
    #     , junction_count

def get_nhd_connections(
    supernetwork_data = {}
    , debuglevel = 0
    , verbose = False
    ):
    if 'mask_file_path' not in supernetwork_data:
        #TODO: doing things this way may mean we are reading the same file twice -- fix this [maybe] by implementing an overloaded return
        supernetwork_data.update({'mask_file_path':None})
        supernetwork_data.update({'mask_layer_string':None})
        supernetwork_data.update({'mask_driver_string':None})
        supernetwork_data.update({'mask_key_col':None})
    return do_connections(
        geo_file_path = supernetwork_data['geo_file_path']
          , key_col = supernetwork_data['key_col']
          , downstream_col = supernetwork_data['downstream_col']
          , length_col = supernetwork_data['length_col']
          , terminal_code = supernetwork_data['terminal_code']
          , title_string = supernetwork_data['title_string']
          , driver_string = supernetwork_data['driver_string']
          , layer_string = supernetwork_data['layer_string']
          , mask_file_path = supernetwork_data['mask_file_path']
          , mask_layer_string = supernetwork_data['mask_layer_string']
          , mask_driver_string = supernetwork_data['mask_driver_string']
          , mask_key_col = supernetwork_data['mask_key_col']
          , debuglevel = debuglevel
          , verbose = verbose
        )

def set_supernetwork_data(
    supernetwork = ''
    , geo_input_folder = None
    , verbose = True
    , debuglevel = 0
    ):

    # The following datasets are extracts from the feature datasets available
    # from https://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/web/data_tools/
    # the CONUS_ge5 and Brazos_LowerColorado_ge5 datasets are included
    # in the github test folder

    supernetwork_options = {
        'Pocono_TEST1'
	,'Pocono_TEST2'
        , 'LowerColorado_Conchos_FULL_RES'
        , 'Brazos_LowerColorado_ge5'
        , 'Brazos_LowerColorado_FULL_RES'
        , 'Brazos_LowerColorado_Named_Streams'
        , 'CONUS_ge5'
        , 'Mainstems_CONUS'
        , 'CONUS_Named_Streams'
        , 'CONUS_FULL_RES_v20'
    }
    if supernetwork not in supernetwork_options:
        print('Note: please call function with supernetworks set to one of the following:')
        for s in supernetwork_options:
            print(f"'{s}'")
        exit()
    elif supernetwork == 'Pocono_TEST1':
        return {
            'geo_file_path' : os.path.join(geo_input_folder
                    , r'PoconoSampleData1' 
                    , r'PoconoSampleRouteLink1.shp')
            , 'key_col' : 18
            , 'downstream_col' : 23
            , 'length_col' : 5
            , 'manningn_col' : 20
            , 'manningncc_col' : 21
            , 'slope_col' : 10
            , 'bottomwidth_col' : 2
            , 'topwidth_col' : 11
            , 'topwidthcc_col' : 12
            , 'MusK_col' : 7
            , 'MusX_col' : 8
            , 'ChSlp_col' : 12
            , 'terminal_code' : 0
            , 'title_string' : 'Pocono Test Example'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }
    elif supernetwork == 'Pocono_TEST2':
        return {
            'geo_file_path' : os.path.join(geo_input_folder
		    , r'PoconoSampleData2' 
                    , r'PoconoRouteLink_testsamp1_nwm_mc.shp')
            , 'key_col' : 18
            , 'downstream_col' : 23
            , 'length_col' : 5
            , 'manningn_col' : 20
            , 'manningncc_col' : 21
            , 'slope_col' : 10
            , 'bottomwidth_col' : 2
            , 'topwidth_col' : 11
            , 'topwidthcc_col' : 12
            , 'MusK_col' : 6
            , 'MusX_col' : 7
            , 'ChSlp_col' : 3
            , 'terminal_code' : 0
            , 'title_string' : 'Pocono Test 2 Example'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'LowerColorado_Conchos_FULL_RES':
        dict = set_supernetwork_data(
                supernetwork = 'CONUS_FULL_RES_v20'
                , geo_input_folder = geo_input_folder
                )
        dict.update({
            'title_string' : 'NHD 2.0 Conchos Basin of the LowerColorado River' #overwrites other title...
              , 'mask_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'masks'
                    , r'LowerColorado_Conchos_FULL_RES.txt')
              , 'mask_driver_string' : r'csv'
              , 'mask_layer_string' : r''
              , 'mask_key_col' : 0
              , 'mask_name_col' : 1 #TODO: Not used yet.
            })
        return dict

    elif supernetwork == 'Brazos_LowerColorado_ge5':
        return {
            'geo_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'NHD_BrazosLowerColorado_Channels.shp')
            , 'key_col' : 2
            , 'downstream_col' : 7
            , 'length_col' : 6
            , 'manningn_col' : 11
            , 'slope_col' : 10
            , 'bottomwidth_col' : 12
            , 'MusK_col' : 7
            , 'MusX_col' : 8
            , 'ChSlp_col' : 13
            , 'terminal_code' : 0
            , 'title_string' : 'NHD Subset including Brazos + Lower Colorado\nNHD stream orders 5 and greater'
            , 'driver_string' : 'ESRI Shapefile'
            , 'layer_string' : 0
          }

    elif supernetwork == 'Brazos_LowerColorado_FULL_RES':
        dict = set_supernetwork_data(
                supernetwork = 'CONUS_FULL_RES_v20'
                , geo_input_folder = geo_input_folder
                )
        dict.update({
            'title_string' : 'NHD 2.0 Brazos and LowerColorado Basins' #overwrites other title...
              , 'mask_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'masks'
                    , r'Brazos_LowerColorado_FULL_RES.txt')
              , 'mask_driver_string' : r'csv'
              , 'mask_layer_string' : r''
              , 'mask_key_col' : 0
              , 'mask_name_col' : 1 #TODO: Not used yet.
            })
        return dict

    elif supernetwork == 'Brazos_LowerColorado_Named_Streams':
        dict = set_supernetwork_data(
                supernetwork = 'CONUS_FULL_RES_v20'
                , geo_input_folder = geo_input_folder
                )
        dict.update({
            'title_string' : 'NHD 2.0 GNIS labeled streams in the Brazos and LowerColorado Basins' #overwrites other title...
              , 'mask_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'masks'
                    , r'Brazos_LowerColorado_Named_Streams.csv')
              , 'mask_driver_string' : r'csv'
              , 'mask_layer_string' : r''
              , 'mask_key_col' : 0
              , 'mask_name_col' : 1 #TODO: Not used yet.
            })
        return dict

    elif supernetwork == 'CONUS_ge5':
        dict = set_supernetwork_data(
                supernetwork = 'CONUS_FULL_RES_v20'
                , geo_input_folder = geo_input_folder
                )
        dict.update({
            'title_string' : 'NHD CONUS Order 5 and Greater' #overwrites other title...
              , 'mask_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'masks'
                    , r'CONUS_ge5.txt')
              , 'mask_driver_string' : r'csv'
              , 'mask_layer_string' : r''
              , 'mask_key_col' : 0
              , 'mask_name_col' : 1 #TODO: Not used yet.
            })
        return dict

    elif supernetwork == 'Mainstems_CONUS':
        return {
            'geo_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'conus_routeLink_subset.nc')
            , 'key_col' : 0
            , 'downstream_col' : 2
            , 'length_col' : 10
            , 'manningn_col' : 11
            , 'manningncc_col' : 20
            , 'slope_col' : 12
            , 'bottomwidth_col' : 14
            , 'topwidth_col' : 22
            , 'topwidthcc_col' : 21
            , 'MusK_col' : 8
            , 'MusX_col' : 9
            , 'ChSlp_col' : 13
            , 'terminal_code' : 0
            , 'title_string' : 'CONUS "Mainstem"'
            , 'driver_string' : 'NetCDF'
            , 'layer_string' : 0
          }

    elif supernetwork == 'CONUS_Named_Streams':
        dict = set_supernetwork_data(
                supernetwork = 'CONUS_FULL_RES_v20'
                , geo_input_folder = geo_input_folder
                )
        dict.update({
            'title_string' : 'CONUS NWM v2.0 only GNIS labeled streams' #overwrites other title...
              , 'mask_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'masks'
                    , r'nwm_reaches_conus_v21_wgnis_name.csv')
              , 'mask_driver_string' : r'csv'
              , 'mask_layer_string' : r''
              , 'mask_key_col' : 0
              , 'mask_name_col' : 1 #TODO: Not used yet.
            })
        return dict

    elif supernetwork == 'CONUS_FULL_RES_v20':
        return {
            'geo_file_path' : os.path.join(geo_input_folder
                    , r'Channels'
                    , r'RouteLink_NWMv2.0_20190517_cheyenne_pull.nc')
            , 'key_col' : 0
            , 'downstream_col' : 2
            , 'length_col' : 10
            , 'manningn_col' : 11
            , 'manningncc_col' : 20
            , 'slope_col' : 12
            , 'bottomwidth_col' : 14
            , 'topwidth_col' : 22
            , 'topwidthcc_col' : 21
            , 'MusK_col' : 8
            , 'MusX_col' : 9
            , 'ChSlp_col' : 13
            , 'terminal_code' : 0
            , 'title_string' : 'CONUS Full Resolution NWM v2.0'
            , 'driver_string' : 'NetCDF'
            , 'layer_string' : 0
          }

def set_networks(
    supernetwork = ''
    , geo_input_folder = None
    , verbose = True
    , debuglevel = 0
    ):

    supernetwork_data = set_supernetwork_data(
      supernetwork = supernetwork
      , geo_input_folder = geo_input_folder
      ) 
    supernetwork_values = get_nhd_connections(
      supernetwork_data = supernetwork_data
      , verbose = verbose
      , debuglevel = debuglevel
      )
    return supernetwork_data, supernetwork_values

