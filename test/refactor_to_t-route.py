#!/usr/bin/env python3

import argparse
import os
import geopandas as gpd
import pandas as pd
import json
import numpy as np
import xarray as xr
from pathlib import Path
from nwm_routing.__main__ import _run_everything_v02
from troute.nhd_io import read_custom_input

'''
Notes for Mike:
    * mising segments in route link file (maybe OCONUS?, adding the 'alt' and 'NHDWaterbodyComID' fields on your end should resolve this)
    * reference hydrofabirc formatting: can we run the same processing script on the reference fabric to generate similar attributes
    * not fully understanding the proportional weighting notaion; I will need this to calculate q_lat time series for the refactored hydrofabric (and currently adding attributes)
TODOs:
    * package scripts into continuous workflow from refactored to t-route
    * create formal template for evaluation metrics
    * add q_lat time series functionality
'''

def generate_qlat(refactored_streams, orig_q_lat):
    # generate q_lat time series for new network from NWM v2.1 IDs

    orig_q_lat = orig_q_lat.rename(columns={orig_q_lat.columns[0]: 'ID'})
    new_qlat =  pd.DataFrame()

    # transpose dataframe
    orig_q_lat_melt = pd.melt(orig_q_lat, id_vars=['ID'], var_name='Date', value_name='q_lat')

    for index, segment in refactored_streams.iterrows():

        crosswalked_segments = segment.comids.split(',')
        q_lat_calc = np.zeros(len(orig_q_lat_melt.Date.unique()))
        new_partial_qlat = pd.DataFrame(columns={'Date','q_lat'}) # columns={'ID', 'Date','q_lat'}

        for s in crosswalked_segments:

            if '.' in s:

                # Get ID and fraction of crosswalk
                partial_id, fraction = s.split('.')

                # Get partial q_lat
                partial_q_lat = orig_q_lat_melt.loc[orig_q_lat_melt.ID==int(partial_id)].copy().reset_index()

                # Calculate partial q_lat
                partial_q_lat["q_lat"] = int(fraction) * partial_q_lat["q_lat"]
                q_lat_calc_temp = list(int(fraction) * partial_q_lat["q_lat"] )
            else:
                # Get partial q_lat
                partial_q_lat = orig_q_lat_melt.loc[orig_q_lat_melt.ID==int(s)].copy().reset_index()
                q_lat_calc_temp = list(partial_q_lat["q_lat"])

            # aggregate q_lat
            q_lat_calc = pd.Series([a + b for a, b in zip(q_lat_calc, q_lat_calc_temp)], name='q_lat')

        new_partial_qlat = new_partial_qlat.append(pd.concat([partial_q_lat.Date, q_lat_calc], axis=1))
        new_partial_qlat.insert(0, 'ID', segment.ID)

        # Add q_lat for new segment ID
        new_qlat = new_qlat.append(new_partial_qlat,ignore_index=True)

    return new_qlat


def recalculate_flows(new_network_flows, q_lat_sample_path):
    # generate Q time series for NWM v2.1 network IDs from troute
    print ('coming soon')


def convert_results_to_df(results, nts, return_courant = False):
    courant = pd.DataFrame()

    qvd_columns = pd.MultiIndex.from_product(
        [range(nts), ["q", "v", "d"]]
    ).to_flat_index()

    flowveldepth = pd.concat(
        [pd.DataFrame(r[1], index=r[0], columns=qvd_columns) for r in results],
        copy=False,
    )

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

    return flowveldepth, courant


def run_troute_from_script(custom_input_file):
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
    ) = read_custom_input(custom_input_file)

    # Modify any parameters from a template file inline here
    run_parameters['nts'] = nts

    # Then run the model and get the results
    results = _run_everything_v02(
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

    return results


def add_routelink_attributes(routlink_csv_filename, stream_layer, routlink_netcdf_filename=None):

    # Convert CONUS routelink NetCDF to CSV
    if not os.path.isfile(routlink_csv_filename):
        routlink_netcdf = xr.open_dataset(routlink_netcdf_filename)
        routlink_df = routlink_netcdf.to_dataframe()
        routlink_df.to_csv(routlink_csv_filename)
    else:
        # Read in routlink parameters
        routlink_df = pd.read_csv(routlink_csv_filename)

    # Add missing attributes
    stream_layer["NHDWaterbodyComID"] = None
    stream_layer["alt"] = None
    fraction_dict = {}
    ids_missing_parameters = []
    comid_set = set()
    for index, segment in stream_layer.iterrows():

        original_segment = segment.comids.split(',')

        # Segments with multiple crosswalk IDs
        if len(original_segment) > 1:

            temp_segs = pd.DataFrame()
            for i in original_segment:

                if '.' in i:
                    # Get ID and fraction of crosswalk
                    i, fraction = i.split('.')

                    if i in fraction_dict.keys():
                        fraction_dict[i] += int(fraction)/10
                    else:
                        fraction_dict[i] = int(fraction)/10

                comid_set.add(i)
                temp_seg = routlink_df.loc[routlink_df.link == int(i)].copy()
                temp_segs = temp_segs.append(temp_seg)

                # Average channel properties
                stream_layer.loc[stream_layer.ID == segment.ID,"NHDWaterbodyComID"] = str(temp_segs['NHDWaterbodyComID'].max())
                stream_layer.loc[stream_layer.ID == segment.ID,"alt"] = np.round(temp_segs['alt'].mean(),2)

        # Segments with single crosswalk ID
        elif len(original_segment) == 1:

            if '.' in original_segment[0]:
                    # Get ID and fraction of crosswalk
                    i, fraction = original_segment[0].split('.')

                    if i in fraction_dict.keys():
                        fraction_dict[i] += int(fraction)/10
                    else:
                        fraction_dict[i] = int(fraction)/10
            else:
                i = original_segment[0]

            comid_set.add(i)
            # Copy channel properties
            temp_seg = routlink_df.loc[routlink_df.link == int(i)]

            if temp_seg.empty:
                print(f"segment ID {int(i)} is not in route link file")
                ids_missing_parameters = ids_missing_parameters + [i]

            else:
                stream_layer.loc[stream_layer.ID == segment.ID,"NHDWaterbodyComID"] = str(temp_seg['NHDWaterbodyComID'].item())
                stream_layer.loc[stream_layer.ID == segment.ID,"alt"] = temp_seg['alt'].item()
        else:
            print (f"no crosswalk for segment {segment.ID}")

    return stream_layer, ids_missing_parameters, comid_set, fraction_dict


def clean_dataset(original_df):

    # Add column
    original_df.index.name = 'COMID'
    original_df = original_df.reset_index()

    # Transpose dataframe
    original_df_melt = pd.melt(original_df, id_vars=['COMID'], var_name='parameters', value_name='values')

    # Unpack parameter column
    original_df_melt[['routing_timestep', 'routing_variable']] = pd.DataFrame(original_df_melt['parameters'].tolist(), index=original_df_melt.index)

    if original_df_melt.routing_timestep.dtype != 'int': original_df_melt.routing_timestep = original_df_melt.routing_timestep.astype(int)

    original_df_melt = original_df_melt.drop(columns=['parameters'])
    original_df_melt = original_df_melt[original_df_melt['routing_variable'] == 'cn']

    return original_df_melt


def replace_line(filename, dictionary, new_filename):

    with open(filename,"r") as f:
        get_all=f.readlines()

    with open(new_filename,'w') as f:
        for i,line in enumerate(get_all,1):
            if ':' in line:
                try:
                    att, var = line.split(':')
                except:
                    att, var, comment = line.split(':')
                att = att.strip()
                if att in yaml_dict.keys():
                    f.writelines(yaml_dict[att])
                else:
                    f.writelines(line)
            else:
                f.writelines(line)


def create_cn_summary_table(courant_results, stream_network, network):

    bad_segments_list = courant_results[courant_results['values']>= 1].COMID.unique()
    violating_segs = stream_network.loc[stream_network.link.isin(bad_segments_list)]

    cn_summary_table = pd.DataFrame({
                                    'Network': [network],
                                    'Min Segment Length (m)': np.round(stream_network.Length.min(),1),
                                    'Max Segment Length (m)': np.round(stream_network.Length.max(),1),
                                    'Avg Segment Length (m)': np.round(stream_network.Length.mean(),1),
                                    'Total Segment Length (m)': np.round(stream_network.Length.sum(),1),
                                    'Total Segments': stream_network.Length.count(),
                                    'Total Timesteps': len(courant_results.routing_timestep.unique()),
                                    'CN Violations': len(courant_results[courant_results['values']>= 1]),
                                    'Number of Segments Violating CN': len(bad_segments_list),
                                    'Length of Network violation CN': violating_segs.Length.sum(),
                                    '% of Network violating CN': np.round(100*(violating_segs.Length.sum()/stream_network.Length.sum()),1)
                                    })

    return cn_summary_table

def run_troute(huc_id,nts,testing_dir):

    '''
    1. Populate a refactored stream network with channel routing properties from NWM v2.1 route link file.
       Subset by available channel routing properties.
    2. Run MC/diffusive wave routing test on the Hurricane Florence Domain to calculate Courant Numbers for each segment for every timestep (separate script)
    3. Clean up CN output format and generate summary statistics

    reporting
    1) CN violations: segment count: timesteps
    2) min, max, average, median stream length for the original and refactored networks
    3) check number of unique segments that violate CN (account for multiple time step violations)
    4) % Length of stream network with CN violations over period of time
    5) network order -- if we could bin by order, that might be interesting
    '''

    testing_dir = Path(testing_dir)
    default_yaml_filename = testing_dir / 'refactored_test.yaml'

    inputs_dir = testing_dir / "inputs"
    routlink_csv_filename = inputs_dir / 'RouteLink_CONUS.csv'
    routlink_netcdf_filename = inputs_dir / 'RouteLink_CONUS.nc'
    nwm_v2_1_streams_filename = inputs_dir / "v2.1" / str("GF_" + str(huc_id) + ".gpkg")
    nwm_v2_1_streams_subset_filename = inputs_dir / "v2.1" / str("routlink_" + str(huc_id) + ".gpkg")
    refactored_streams_dir = inputs_dir / "refactored"
    refactored_hydrofabrics = os.listdir(refactored_streams_dir)

    outputs_dir = testing_dir / "outputs"
    if not Path(outputs_dir).is_dir():
        os.mkdir(Path(outputs_dir))
    aggregate_cn_table_filename = outputs_dir / 'aggregate_cn_table.csv'
    aggregate_cn_summary_table_filename = outputs_dir / 'aggregate_cn_summary_table.csv'

    # Calculate base metrics
    if not os.path.isfile(nwm_v2_1_streams_subset_filename):

        nwm_v2_1_streams = gpd.read_file(nwm_v2_1_streams_filename,layername='nhd_flowline')
        subset_ids = nwm_v2_1_streams.COMID.to_list()

        routlink_df = pd.read_csv(routlink_csv_filename)

        subset_routelink = routlink_df.loc[routlink_df.link.isin(subset_ids)].copy()
        subset_routelink.reset_index(drop=True, inplace=True)
        subset_routelink.drop(columns=['feature_id'], inplace=True)

        subset_routelink = nwm_v2_1_streams[['COMID','geometry']].merge(subset_routelink, left_on=['COMID'], right_on=['link'])

        subset_routelink.to_file(nwm_v2_1_streams_subset_filename,layer='nwm_flowline',driver='GPKG')

        ## Create t-route input file
        routlink_cn_csv_filename = os.path.join(testing_dir,'courant_routlink_' + huc_id + '.csv')
        output_filename = os.path.join(os.path.dirname(default_yaml_filename), 'reference_' + huc_id + '.yaml')

        yaml_dict = {}
        yaml_dict["title_string"] = "    title_string: " +  "routlink_" + huc_id + "\n"
        yaml_dict["geo_file_path"] = '    geo_file_path: ' +  '\"' + str(nwm_v2_1_streams_subset_filename) + '\"' + '\n'
        yaml_dict["layer_string"] = '    layer_string: nhd_flowline' + '\n'

        # Write yaml file
        replace_line(default_yaml_filename, yaml_dict, output_filename)

        ## Run MC/diffusive wave routing test to calculate Courant Numbers for each segment for every timestep (calls compiled Fortran module)
        # Run t-route
        t_route_results = run_troute_from_script(str(output_filename))
        # Convert results to dataframe
        fvd, courant = convert_results_to_df(t_route_results, nts, True)

        # Clean up Courant Number output table
        tidy_network_cn = clean_dataset(courant)
        tidy_network_cn['parameters'] = 'nwm_flowline'

        network = 'reference_' + huc_id
        cn_summary_table = create_cn_summary_table(tidy_network_cn, subset_routelink, network)

        cn_summary_table.to_csv(aggregate_cn_summary_table_filename,index=False)

        del nwm_v2_1_streams, subset_routelink, t_route_results, fvd, courant, tidy_network_cn

        # # Write out CN values to aggregate table
        # if os.path.isfile(aggregate_cn_table_filename):
        #     tidy_network_cn.to_csv(aggregate_cn_table_filename,index=False, mode='a',header=False)
        # else:
        #     tidy_network_cn.to_csv(aggregate_cn_table_filename,index=False)

    # Calculate metrics for each refactored hydrofabrics
    for hydrofabric in refactored_hydrofabrics:
        print (f"running hydrofabric {hydrofabric}")

        # Read in refactored stream network
        refactored_streams = gpd.read_file(Path(refactored_streams_dir) / hydrofabric)

        if not set(["NHDWaterbodyComID","alt"]).issubset(refactored_streams.columns):
            # Populate stream network with channel routing properties from NWM v2.1 route link file.
            refactored_streams, ids_missing_parameters, comid_set, fraction_dict = add_routelink_attributes(routlink_csv_filename, refactored_streams, routlink_netcdf_filename)

            with open(Path(testing_dir / str(hydrofabric + '_comid_set.lst')), 'w') as f:
                for comid in list(comid_set):
                    f.write("%s\n" % comid)

            with open(Path(testing_dir / str(hydrofabric + '_comid_set.json')), 'w') as f:
                json.dump(fraction_dict, f)

            print(f"{len(ids_missing_parameters)} segments were not in routelink file")

            # Remove NA segments and clean up column names (this can be improved)
            refactored_streams = refactored_streams.dropna(subset=['n','toID'])
            refactored_streams = refactored_streams.rename(columns={'lengthm': 'Length', 'ID': 'link', 'toID': 'to'})

            if refactored_streams.to.dtype != 'int': refactored_streams.to = refactored_streams.to.astype(int)
            if refactored_streams.NHDWaterbodyComID.dtype != 'int': refactored_streams.NHDWaterbodyComID = refactored_streams.NHDWaterbodyComID.astype(int)

            # Save updated gpkg
            clean_up_dir = Path(os.path.dirname(refactored_streams_dir)) / "refactored_cleaned_up"
            if not clean_up_dir.is_dir():
                os.mkdir(clean_up_dir)
            refactored_streams_filename = Path(clean_up_dir) / hydrofabric
            refactored_streams.to_file(refactored_streams_filename, layer='refactored',driver='GPKG')

        else:
            refactored_streams_filename = Path(refactored_streams_dir) / hydrofabric

        ## Create t-route input file
        yaml_dict = {}
        hydrofabric_version, layer_type = hydrofabric.split('.')
        yaml_dict["title_string"] = "    title_string:" +  " refactored_" + hydrofabric_version + "\n"
        yaml_dict["geo_file_path"] = '    geo_file_path: ' +  '\"' + str(refactored_streams_filename) + '\"' + '\n'
        yaml_dict["layer_string"] = '    layer_string: refactored' + '\n'
        # yaml_dict["qlat_const:"] = "0.1"
        # yaml_dict["qlat_input_file:"] = "../../test/input/florence_933020089/FORCING/SeptemberFlorenceQlateral_wDates.csv"

        output_filename = os.path.join(os.path.dirname(default_yaml_filename), 'refactored_' + hydrofabric_version + '.yaml')

        # Write yaml file
        replace_line(default_yaml_filename, yaml_dict, output_filename)

        ## Run MC/diffusive wave routing test to calculate Courant Numbers for each segment for every timestep (calls compiled Fortran module)
        # Run t-route
        t_route_results = run_troute_from_script(str(output_filename))
        # Convert results to dataframe
        fvd, courant = convert_results_to_df(t_route_results, nts, True)

        ## Clean up Courant Number output format and generate summary statistics
        tidy_network_cn = clean_dataset(courant)
        tidy_network_cn['parameters'] = hydrofabric_version

        # rpu, split_flines_meters, collapse, split_main_collapse = hydrofabric_version.split('_')
        # tidy_network_cn['rpu'] = rpu
        # tidy_network_cn['split_flines_meters'] = split_flines_meters
        # tidy_network_cn['collapse'] = collapse
        # tidy_network_cn['split_main_collapse'] = split_main_collapse

        # # Write out CN values to aggregate table
        # if os.path.isfile(aggregate_cn_table_filename):
        #     tidy_network_cn.to_csv(aggregate_cn_table_filename,index=False, mode='a')
        # else:
        #     tidy_network_cn.to_csv(aggregate_cn_table_filename,index=False)

        cn_summary_table = create_cn_summary_table(tidy_network_cn, refactored_streams, hydrofabric_version)

        if os.path.isfile(aggregate_cn_summary_table_filename):
            cn_summary_table.to_csv(aggregate_cn_summary_table_filename,index=False, mode='a',header=False)
        else:
            cn_summary_table.to_csv(aggregate_cn_summary_table_filename,index=False)

        del refactored_streams, t_route_results, fvd, courant, tidy_network_cn

if __name__ == '__main__':

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Run t-route routing algorithms and calculate metrics')
    parser.add_argument('-huc', '--huc-id', help='HUC2 ID',required=True,type=str)
    parser.add_argument('-nts', '--nts', help='number of timesteps',required=True,type=int)
    parser.add_argument('-dir', '--testing-dir', help='testing directory',required=True)

    # Extract to dictionary and assign to variables.
    args = vars(parser.parse_args())

    run_troute(**args)
