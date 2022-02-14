reservoir_ids = [401, 402, 403]
network_clean = [
    [0, 456, -999, 0],
    [1, 178, 4, 0],
    [2, 394, 0, 0],
    [3, 301, 2, 0],
    [4, 798, 0, 403],
    [5, 679, 4, 403],
    [6, 523, 0, 0],
    [7, 815, 2, 0],
    [8, 841, -999, 0],
    [9, 514, 8, 0],
    [10, 458, 9, 0],
    [11, 832, 10, 0],
    [12, 543, 11, 0],
    [13, 240, 12, 0],
    [14, 548, 13, 0],
    [15, 920, 14, 0],
    [16, 920, 15, 401],
    [17, 514, 16, 401],
    [18, 458, 17, 0],
    [180, 458, 17, 0],
    [181, 458, 180, 0],
    [19, 832, 18, 0],
    [20, 543, 19, 0],
    [21, 240, 16, 401],
    [22, 548, 21, 0],
    [23, 920, 22, 0],
    [24, 240, 23, 0],
    [25, 548, 12, 0],
    [26, 920, 25, 402],
    [27, 920, 26, 402],
    [28, 920, 27, 0],
    [2800, 920, 2700, 0],
]
# Should create independent networks with terminal nodes at 2800, 0, and 8

# 50-21, 60-62, 70-73, 81-84 are sets of circular networks, 
# which is an error and should be filtered
network_circulars = [
    [50, 178, 51, 0],
    [51, 178, 50, 0],
    [60, 178, 61, 0],
    [61, 178, 62, 0],
    [62, 178, 60, 0],
    [70, 178, 71, 0],
    [71, 178, 72, 0],
    [72, 178, 73, 0],
    [73, 178, 70, 0],
    [80, 178, 81, 0],
    [81, 178, 82, 0],
    [82, 178, 83, 0],
    [83, 178, 84, 0],
    [84, 178, 80, 0],
    [0, 456, -999, 0],
    [1, 178, 4, 0],
    [2, 394, 0, 0],
    [3, 301, 2, 0],
    [4, 798, 0, 0],
    [5, 679, 4, 0],
    [6, 523, 0, 0],
    [7, 815, 2, 0],
    [8, 841, -999, 0],
    [9, 514, 8, 0],
    [10, 458, 9, 0],
    [11, 832, 10, 0],
    [12, 543, 11, 0],
    [13, 240, 12, 0],
    [14, 548, 13, 0],
    [15, 920, 14, 0],
    [16, 920, 15, 401],
    [17, 514, 16, 401],
    [18, 458, 17, 0],
    [180, 458, 17, 0],
    [181, 458, 180, 0],
    [19, 832, 18, 0],
    [20, 543, 19, 0],
    [21, 240, 16, 401],
    [22, 548, 21, 0],
    [23, 920, 22, 0],
    [24, 240, 23, 0],
    [25, 548, 12, 0],
    [26, 920, 25, 0],
    [27, 920, 26, 0],
    [28, 920, 27, 0],
    [2800, 920, 2700, 0], 
]

test_key_col = 0
test_downstream_col = 2
test_length_col = 1
test_terminal_code = -999
test_waterbody_col = 3
test_waterbody_null_code = 0

test_columns = {
    "key": 0,
    "dx": 1,
    "downstream": 2,
    "waterbody": 3,
}

reverse_test_columns = {0: "key", 1: "dx", 2: "downstream", 3: "waterbody"}

expected_connections = {
    0: [],
    1: [4],
    2: [0],
    3: [2],
    4: [0],
    5: [4],
    6: [0],
    7: [2],
    8: [],
    9: [8],
    10: [9],
    11: [10],
    12: [11],
    13: [12],
    14: [13],
    15: [14],
    16: [15],
    17: [16],
    18: [17],
    180: [17],
    181: [180],
    19: [18],
    20: [19],
    21: [16],
    22: [21],
    23: [22],
    24: [23],
    25: [12],
    26: [25],
    27: [26],
    28: [27],
    2800: [],
}

expected_rconn = {0: [2, 4, 6], 1: [], 4: [1, 5], 2: [3, 7], 3: [], 5: [], 6: [], 7: [], 8: [9], 9: [10], 10: [11], 11: [12], 12: [13, 25], 13: [14], 14: [15], 15: [16], 16: [17, 21], 17: [18, 180], 18: [19], 180: [181], 181: [], 19: [20], 20: [], 21: [22], 22: [23], 23: [24], 24: [], 25: [26], 26: [27], 27: [28], 28: [], 2800: []}

expected_wbody_connections = {4: 403, 5: 403, 16: 401, 17: 401, 21: 401, 26: 402, 27: 402}

# network_with_substituted_reservoirs = []
# rconn = {"down":["up (is a reservoir)"]}
# path_func = partial(nhd_network.split_at_waterbodies_and_junctions, ["up (is a reservoir)"], rconn)
# nhd_network.dfs_decomposition_depth_tuple(rconn, path_func)
# list(nhd_network.dfs_decomposition_depth_tuple(rconn, path_func))
# [(0, ['down'])]

import pandas as pd
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_network as nhd_network

test_param_df = pd.DataFrame(network_clean)
test_param_df = test_param_df.rename(columns=nhd_network.reverse_dict(test_columns))
test_param_df = test_param_df.set_index("key")


def test_reverse_dict():

    result = nhd_network.reverse_dict(test_columns)
    assert result == reverse_test_columns


def test_build_connections():

    # There can be an externally determined terminal code -- that's this first value
    terminal_codes = set()
    terminal_codes.add(test_terminal_code)
    # ... but there may also be off-domain nodes that are not explicitly identified
    # but which are terminal (i.e., off-domain) as a result of a mask or some other
    # an interior domain truncation that results in a
    # otherwise valid node value being pointed to, but which is masked out or
    # being intentionally separated into another domain.
    terminal_codes = terminal_codes | set(
        test_param_df[~test_param_df["downstream"].isin(test_param_df.index)][
            "downstream"
        ].values
    )

    connections = nhd_network.extract_connections(
        test_param_df, "downstream", terminal_codes
    )
    assert connections == expected_connections


def test_reverse_network():
    connections = expected_connections
    rconn = nhd_network.reverse_network(connections)
    assert expected_rconn == rconn
    rrconn = nhd_network.reverse_network(rconn)
    assert rrconn == connections


def test_extract_waterbodies():
    wbody_connections = nhd_network.extract_waterbody_connections(test_param_df, "waterbody", test_waterbody_null_code)
    assert wbody_connections == expected_wbody_connections

