import troute.nhd_io as nhd_io 

def ngen_preprocess():
    """

    """

    ngen_network_df = nhd_io.read_geopandas( "flowpath_data.geojson" )
    print(ngen_network_df) 
    # Extract subset here if needed

    # map the connections into dictionary form
    ngen_network_dict = dict(zip(ngen_network_df.id, ngen_network_df.toid))
    
    def node_key_func(x):
        return int(x[3:])

    # Extract the ID integer values
    connections = {node_key_func(k): node_key_func(v) for k, v in ngen_network_dict.items()}

    # Convert dictionary connections to data frame and make ID column the index
    connection_df = pd.DataFrame.from_dict(connections, orient='index', columns=['to'])
    # Sort ID index column
    connection_df = connection_df.sort_index()

    connection_df = nhd_io.replace_downstreams(waterbody_df, "to", 0)

    connections = nhd_network.extract_connections(waterbody_df, "to")


    # Channel parameters
    param_df = nhd_id.read_netcdf("RouteLink_NHDPLUS.nc")
    param_df['dt'] = 300.0

    param_df.set_index("link", inplace=True)

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

    routelink_cols = dict([(value,key) for key, value in routelink_cols.items() ])

    param_df.rename(columns=routelink_cols, inplace=True)

