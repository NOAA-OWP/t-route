import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io 
import pandas as pd
import json


# Each catchment has one or more segments/comids from the crosswalk-mapping. 
# Each cacthment has one downstream nexus from the flowpath_data.geojson.
# From the Route Link file, each segment/comid has one downstream segment/comid. 
# From filtering, each catchment can then reduce to having one downstream segment/comid.
# We can then map the single nexus to one downstream segment/comid.
# In the end, multiple nexuses can point to a single downstream segment/comid,
# but a single nexus cannot have multiple downstream segments/comids.

# Actually a single nexus can have mult ds comids. Need to revisit this later


def ngen_preprocess():
    """

    """

    #ngen_network_df = nhd_io.read_geopandas( "flowpath_data.geojson" )
    ngen_network_df = nhd_io.read_geopandas( "/apd_common/anthro/david_code/t-route/test/input/next_gen/flowpath_data.geojson" )
    #ngen_network_df = nhd_io.read_geopandas( "/apd_common/anthro/david_code/ngen/data/flowpath_data.geojson" )
    print("ngen_network_df1") 
    print(ngen_network_df) 
    # Extract subset here if needed

    # map the connections into dictionary form
    #ngen_network_dict = dict(zip(ngen_network_df.id, ngen_network_df.toid))
    ngen_network_dict = dict(zip(ngen_network_df.id, ngen_network_df.to))
    
    def node_key_func(x):
        return int(x[4:])

    # Extract the ID integer values
    connections = {node_key_func(k): node_key_func(v) for k, v in ngen_network_dict.items()}
    print ('connections1')
    print (connections)

    # Convert dictionary connections to data frame and make ID column the index
    connection_df = pd.DataFrame.from_dict(connections, orient='index', columns=['to'])
    print ("connection_df1")
    print (connection_df)
    # Sort ID index column
    connection_df = connection_df.sort_index()
    print ("connection_df2")
    print (connection_df)

    #No longer need this because not setting up catchment network for t-route
    #connection_df = nhd_io.replace_downstreams(connection_df, "to", 0)

    #Set index name to cat-id
    connection_df.index.names = ["cat-id"]

    print ("connection_df3")
    print (connection_df)


    connections = nhd_network.extract_connections(connection_df, "to")

    print ("connections")
    print (connections)


    # Channel parameters
    #param_df = nhd_io.read_netcdf("RouteLink_NHDPLUS.nc")
    #param_df = nhd_io.read_netcdf("............/test_persistence/v3.0_enhancements_model_runs/nwm_v3/model_inputs/NWM/DOMAIN/RouteLink_NWMv2.1.nc")
    param_df = nhd_io.read_netcdf("/apd_common/anthro/backup/water_management_v3.0_archive/nwm_v3/DOMAIN/RouteLink_NWMv2.1.nc")
    param_df['dt'] = 300.0

    param_df.set_index("link", inplace=True)

    #Change index name from link to comid in order to join with the connection_df on its multi-index   
    param_df.index.names = ["comid"]

    #Convert to names used by t-route
    routelink_cols = {
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
        "cs": "ChSlp",
    }

    #Convert to names used by t-route
    routelink_cols = dict([(value,key) for key, value in routelink_cols.items() ])

    #Convert to names used by t-route
    param_df.rename(columns=routelink_cols, inplace=True)

    print ("param_df")
    print (param_df)

    #with open(next_gen_input_folder/'coarse/crosswalk.json') as f:
    #with open(".........../ngen_pybind/crosswalk_mapping/crosswalk-mapping.json") as f:
    with open("/apd_common/anthro/david_code/t-route/test/input/next_gen/crosswalk-mapping.json") as f:
        crosswalk_data = json.load(f)
    #connection_df['comid'] = connection_df.apply(lambda x: crosswalk_data['cat-' + str(x.name)]['outlet_COMID'], axis=1)
    connection_df['comid'] = connection_df.apply(lambda x: crosswalk_data['cat-' + str(x.name)]['COMID'], axis=1)
   
    print ("crosswalk_data")
    #print (crosswalk_data)

    print (list(crosswalk_data.items())) 
    print ("end crosswalk_data -------------------------------")

    cat_id_and_comid_tuple_list = []
    comid_list = []

    for cat_id, value1 in crosswalk_data.items():
        for key2, value2 in value1.items():
            #print (value2)
            for comid in value2:
                #print (comid)
                #cat_id_and_comid_tuple = (cat_id[4:], comid)
                cat_id_and_comid_tuple = (int(cat_id[4:]), int(comid))
                #print (cat_id_and_comid_tuple)
                cat_id_and_comid_tuple_list.append(cat_id_and_comid_tuple)
                
                #don't add duplicates
                if comid not in comid_list:
                    comid_list.append(comid)

    print ("cat_id_and_comid_tuple_list")
    print (cat_id_and_comid_tuple_list)

    print ("comid_list")
    print (comid_list)

    comid_df = pd.DataFrame(comid_list)

    print ("comid_df")
    print (comid_df)

    # Create mask as input to filter sub-domain for t-route
    comid_df.to_csv('ngen_comid_mask.csv', header=False, index=False)


    mult_index = pd.MultiIndex.from_tuples(cat_id_and_comid_tuple_list, names=["cat-id", "comid"])

    print ("mult_index")
    print (mult_index)

    cat_id_and_comid_tuple_df = pd.DataFrame(comid_list, index=mult_index)

    print ("cat_id_and_comid_tuple_df.dtypes")
    print (cat_id_and_comid_tuple_df.dtypes)
    print ("connection_df.dtypes")
    print (connection_df.dtypes)

    print ("cat_id_and_comid_tuple_df")
    print (cat_id_and_comid_tuple_df)

    print ("------------------")

    #Probably delete
    ################################# 
    cat_id_index_df = cat_id_and_comid_tuple_df.index.get_level_values("cat-id")
    print ("======================")
    print ("cat_id_index_df")
    print (cat_id_index_df)

    #print ("cat_id_index_df.loc[62031]")
    #print (cat_id_index_df.loc[62031].values)

    print ("connection_df.loc[cat_id_index_df].values")
    print (connection_df.loc[cat_id_index_df].values)
    #connection_df_with_segs = connection_df.loc[cat_id_index_df].values
    #connection_df = connection_df.set_index(["cat-id"]) 
    #comid_index_df = cat_id_and_comid_tuple_df.index.get_level_values("comid")
    #print ("param_df.loc[10297950]")
    #print (param_df.loc[10297950])
    #cat_id_and_comid_tuple_df["to"] = connection_df.loc[cat_id_index_df].values
    #print ("======================")
    #print ("cat_id_and_comid_tuple_df") 
    #print (cat_id_and_comid_tuple_df) 
    #waterbody_df = waterbody_df.join(nhd_routelink, on='comid', how='left')
    #connection_df = connection_df.join(param_df, on='comid', how='left')
    #print (connection_df.loc[61565])
    #print (connection_df.loc[61875])
    ##################################


    #Join connection_df with the multi-index df 
    connection_df_with_segs = cat_id_and_comid_tuple_df.join(connection_df, how='inner')

    print ("connection_df_with_segs")
    print (connection_df_with_segs)

    #Join the cat-id and seg id df with the param_df
    connection_df_with_segs_and_params = connection_df_with_segs.join(param_df, how='inner')

    print ("connection_df_with_segs_and_params")
    print (connection_df_with_segs_and_params)

    del param_df


    #connection_and_params_df_indexed_by_nexus = connection_df_with_segs_and_params


    #connection_df_with_segs_and_params.to_json("routing_input_parms_to_ngen.json", orient='split', index=False)


    connection_df_with_segs_and_params = connection_df_with_segs_and_params.drop(columns='time')
    connection_df_with_segs_and_params = connection_df_with_segs_and_params.drop(columns='gages')



    connection_dict_with_segs_and_params = connection_df_with_segs_and_params.groupby(level=0).apply(lambda connection_df_with_segs_and_params: connection_df_with_segs_and_params.xs(connection_df_with_segs_and_params.name).to_dict()).to_dict()


    #downstream_comid_to_nexus_id_dict = {}
    nexus_id_to_downstream_comid_dict = {}

    for cat_id_key, params_value in connection_dict_with_segs_and_params.items():
        for comid_key, nexus_id_value in params_value['to'].items():
            print (str(comid_key) + " " + str(nexus_id_value))
            nexus_id = nexus_id_value

        comid_key_list = []
        downstream_value_list = []

        for comid_key, downstream_value in params_value['downstream'].items():
           comid_key_list.append(comid_key)
           downstream_value_list.append(downstream_value)

        downstream_comid = None

        for downstream_value in downstream_value_list:
           if downstream_value not in comid_key_list:
               downstream_comid = downstream_value

        if downstream_comid == None:
           raise ValueError('Did not find a downstream comid for catchment: ' + str(cat_id_key))

        print ("Nexus " + str(nexus_id) + ": DS COMID " + str(downstream_comid))

        #Add key/value of nexus and ds comid to a dict. What happens when a duplicate key/value is entered
        #due to multiple nexuses from mult cats?? Maybe reduce to a set of single key nexus. Or could just check
        #if the key/value is in the dict already and if so then don't add. Consider throwing an error
        #if they key is in the dict but the value for that key in the dict does not match the value
        #that is considering to be added


        #Original thought of dict with ds_comids as keys and nex_is as vals
        #if downstream_comid in downstream_comid_to_nexus_id_dict.keys():
        #    if downstream_comid_to_nexus_id_dict[downstream_comid] != nexus_id:
        #        raise ValueError("COMID: " + str(downstream_comid) + 
        #        " is assigned as downstream to multiple nexus ids")
        print ("nexus map: " + str(nexus_id) + " " + str(downstream_comid))    
 
        if nexus_id in nexus_id_to_downstream_comid_dict.keys():
            if nexus_id_to_downstream_comid_dict[nexus_id] != downstream_comid:
                #raise ValueError("Nexus_id: " + str(nexus_id) + 
                #" is assigned to have multiple downstream segments/comids")
               
                print ("Nexus_id: " + str(nexus_id) +
                " is assigned to have multiple downstream segments/comids")


 
        #Dict of ints right now for both ids
        #downstream_comid_to_nexus_id_dict[downstream_comid] = nexus_id
        nexus_id_to_downstream_comid_dict[nexus_id] = downstream_comid



    print ("connection_dict_with_segs_and_params")
    print (connection_dict_with_segs_and_params)

    #output_file = "routing_input_parms_to_ngen.json"
    #output_file = "downstream_comid_to_nexus_id_mapping.json"
    output_file = "nexus_id_to_downstream_comid_mapping.json"

    #with open(output_file, "w") as open_output_file: 
    #  json.dump(connection_dict_with_segs_and_params, open_output_file, indent=4)

    with open(output_file, "w") as open_output_file: 
      json.dump(nexus_id_to_downstream_comid_dict, open_output_file, indent=4)


if __name__ == "__main__":
    #main()
    ngen_preprocess()

