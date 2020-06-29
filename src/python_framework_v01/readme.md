
### nhd_network_utilities Glossary of objects
and
### nhd_reach_utilities Glossary of objects

##### supernetwork_data
* Basic information about the dataset provided to create the network. 
* This dataset is the Route-Link file for the National Water Model. 
* The 'X_col' attributes in the dictionary indicate which columns of the dataset contain data needed for the routing computation. E.g., the `slope_col` is the column index for the table column with bed slope values. 
* Used by the function `get_nhd_connections`, which calls `do_connections` with the various keys passed to allow development of the connection object and key lists.
e.g., 
```
{'CONUS_ge5':{'geo_file_path': '***REMOVED***user/git/t-route/src/test/input/geo/Channels/NHD_Conus_Channels.shp'
              , 'key_col': 1
              , 'downstream_col': 6
              , 'length_col': 5
              , 'manningn_col': 10
              , 'slope_col': 9
              , 'bottomwidth_col': 11
              , 'MusK_col': 7
              , 'MusX_col': 8
              , 'ChSlp_col': 12
              , 'terminal_code': 0
              , 'title_string': 'NHD CONUS Order 5 and Greater'
              , 'driver_string': 'ESRI Shapefile'
              , 'layer_string': 0}}
```
##### supernetwork_values
the return from the `build_connections_object` function called by `do_connections`
Returns
    #  connections, all_keys, ref_keys, headwater_keys \
    #     , terminal_keys, terminal_ref_keys \
    #     , circular_keys, junction_keys \
    #     , visited_keys, visited_terminal_keys \
    #     , junction_count

###### connections
Each connection is keyed with an nwm/nhd segment id and contains references to the downstream and upstream segments, the length, and the full data row from the original table. 
{328271: {'downstream': 328275
           , 'length': 6685.0
           , 'data': [2283, 328271, 388099, 328271, 8, 6685.0
                      , 328275, 3600.0, 0.2, 0.0, 0.045, 85.6765
                      , 0.05, 0.0, None, -9999, 27.1069, -99.4416, 92.41, 0
                      , 1337104, 0.0621483820102, <shapely.geometry.linestring.LineString object at 0x7fec974f7a58>]
           , 'upstreams': {1131002455}}}
'upstreams' for headwaters are a single valued with the terminal code i.e., for a headwater segment 20337544, 
```
>>> print(connections[20337544]['upstreams'])
```
yields 
```
{0}
```

#####networks
The `compose_networks` function takes the independent network outlets identified while building the connections object, and analyses each independent network to break it into reaches. A reach is defined by any two of a headwater, a junction, or a terminal segment (the outlet into the ocean or into an interior sink). The reaches are assigned a raw sequential order as a distance from the terminal segment and the upstream reaches and downstream reach are added as keyed values.
<br><br>
The example below is a one-reach, 50-segment network that appears in the CONUS_ge5 dataset. The `'junctions'` element of the dictionary is an empty `set()` becuase there are no bifurcations. The set would normally contain the reaches with more than one upstream reach -- which is any reach that is not a headwater by the definition we have used.
```
{20331522: {'total_segment_count': 50
            , 'total_junction_count': 0
            , 'maximum_reach_order': 0
            , 'junctions': set()
            , 'headwater_reaches': {20337544}
            , 'reaches': {20337544: {'reach_tail': 20331522
                                      , 'downstream_reach': 0
                                      , 'reach_head': 20337544
                                      , 'order': 0
                                      , 'segments': {20331520, 20334592, 20331522, 20337538, 20331398, 20337544, 20331402, 948070410, 20337546, 948070157, 948070158, 20332816, 20332818, 948070163, 948070164, 948070165, 948070166, 948070167, 948070168, 948070169, 948070170, 948070171, 20332824, 948070173, 20331934, 20332828, 948070176, 948070177, 948070394, 948070178, 948070179, 20332832, 20331434, 948070396, 20332858, 20332860, 20331388, 20331980, 20331982, 20331984, 20331990, 20332120, 20331486, 20331490, 20331498, 20331500, 20331512, 948070393, 20331514, 20331516}
                                      , 'upstream_reaches': {0}}}
            , 'terminal_reach': 20337544}}
```

