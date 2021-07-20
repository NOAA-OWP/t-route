import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io 
import pandas as pd
import json


def ngen_preprocess():
    """

    """

    #ngen_network_df = nhd_io.read_geopandas( "flowpath_data.geojson" )
    ngen_network_df = nhd_io.read_geopandas( "/glade/work/dmattern/ngen_pybind/crosswalk_mapping/flowpath_data.geojson" )
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

    connection_df = nhd_io.replace_downstreams(connection_df, "to", 0)

    #Set index name to cat-id
    connection_df.index.names = ["cat-id"]

    print ("connection_df3")
    print (connection_df)


    connections = nhd_network.extract_connections(connection_df, "to")

    print ("connections")
    print (connections)


    # Channel parameters
    #param_df = nhd_io.read_netcdf("RouteLink_NHDPLUS.nc")
    param_df = nhd_io.read_netcdf("/glade/scratch/mehdi/test_persistence/v3.0_enhancements_model_runs/nwm_v3/model_inputs/NWM/DOMAIN/RouteLink_NWMv2.1.nc")
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
    with open("/glade/work/dmattern/ngen_pybind/crosswalk_mapping/crosswalk-mapping.json") as f:
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
                comid_list.append(comid)

    print ("cat_id_and_comid_tuple_list")
    print (cat_id_and_comid_tuple_list)

    print ("comid_list")
    print (comid_list)

    comid_df = pd.DataFrame(comid_list)

    print ("comid_df")
    print (comid_df)

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


if __name__ == "__main__":
    #main()
    ngen_preprocess()

