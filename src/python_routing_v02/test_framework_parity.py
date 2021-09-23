from nwm_routing.__main__ import _run_everything_v02  # , _input_handler_v02
# from nwm_routing.__main__ import _run_everything_v03  # , _input_handler_v03
# from troute.routing.compute import compute_nhd_routing_v02
from troute.nhd_io import read_custom_input
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
# import subprocess
# import os
# result = os.system("python -m nwm_routing -V2 -f ../../test/input/yaml/CustomInput.yaml")
#
# version = "-V2"
# subp_c = ["python", "-m", "nwm_routing", version, "-f", test_file]
# result = subprocess.call(subp_c)
# 
# import pickle
# with open("test_data.pickle","rb") as f:
#     input_tuple = pickle.load(f)
# results = compute_nhd_routing_v02(*input_tuple)


def _produce_result_set(results):
    num_networks = len(results)
    id_list_j = 0
    values_list_j = 1
    file_base = "custom_test_standard_result"
    test_folder = Path("../../test/output/Pocono_TEST1/standard_test")

    for net_i in range(num_networks):
        ext_str_id = f".net{net_i:05d}_id.out" 
        ext_str_values = f".net{net_i:05d}_values.out" 
        with open(test_folder / (file_base + ext_str_id), "wb") as fid:
            np.save(fid, results[net_i][id_list_j], allow_pickle = False)
        with open(test_folder / (file_base + ext_str_values), "wb") as fid:
            np.save(fid, results[net_i][values_list_j], allow_pickle = False)


def test_run_everything_v02():
    """
    Integration test for python_routing_v02. List of v02 YAML files with
    final flow of the parity node.
    """

    # The hybrid usgs AnA below does not run correctly whenever the
    # above hybrid usgs is run first. This is because the Fortran
    # singleton timeslice reader is initialized in memory from the 
    # above call with an input timeslice directory set, and that
    # timeslice reader persists in memory for the below call. The
    # below call when run by itself initializes the timeslice reader
    # to a different input timeslice directory, but when run after
    # the above call, the timeslice reader will still search for 
    # timeslice files in the directory specified in the above call.
    input_file_and_result_tuple_list = [
        ("../../test/input/yaml/Croton_NY_levelpool.yaml", 1.415740),
        ("../../test/input/yaml/Croton_NY_hybrid_usgs.yaml", 1.084535),
        #("../../test/input/yaml/Croton_NY_hybrid_usgs_AnA.yaml", 1.3),
        ("../../test/input/yaml/Croton_NY_rfc.yaml", 8.681314),
        # TODO: Updates needed to verify parity for Croton_NY_levelpool_diffusive.yaml
        #("../../test/input/yaml/Croton_NY_levelpool_diffusive.yaml", 0.0)
        ("../../test/input/yaml/CustomInput.yaml", 0.870077),
        ("../../test/input/yaml/Florence_Benchmark.yaml", 0.0043528992)
        #("../../test/input/yaml/Florence_Benchmark_da.yaml", 0.0)
    ]

    for input_file_and_result_tuple in input_file_and_result_tuple_list:

        custom_input_file = input_file_and_result_tuple[0]
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

        nts = run_parameters["nts"] 
        dt = run_parameters["dt"] 

        compare_node = parity_parameters["parity_check_compare_node"]

        # construct a dataframe of simulated flows
        fdv_columns = pd.MultiIndex.from_product([range(nts), ["q", "v", "d"]])
        flowveldepth = pd.concat(
            [pd.DataFrame(r[1], index=r[0], columns=fdv_columns) for r in results],
            copy=False,
        )
        flowveldepth = flowveldepth.sort_index()

        flows = flowveldepth.loc[:, (slice(None), "q")]
        flows = flows.T.reset_index(level=[0, 1])
        flows.rename(columns={"level_0": "Timestep", "level_1": "Parameter"}, inplace=True)
        flows["Time (d)"] = ((flows.Timestep + 1) * dt) / (24 * 60 * 60)
        flows = flows.set_index("Time (d)")

        depths = flowveldepth.loc[:, (slice(None), "d")]
        depths = depths.T.reset_index(level=[0, 1])
        depths.rename(columns={"level_0": "Timestep", "level_1": "Parameter"}, inplace=True)
        depths["Time (d)"] = ((depths.Timestep + 1) * dt) / (24 * 60 * 60)
        depths = depths.set_index("Time (d)")

        # Construct dataframe for results from set compare_node.
        # TODO: Consider testing reservoir water elevation also 
        trt = pd.DataFrame(
            flows.loc[:, compare_node].values,
            columns=["flow, t-route (cms)"],
        )

        final_flow = trt.iloc[-1][0]

        expected_final_flow = input_file_and_result_tuple[1]
        
        assert expected_final_flow == pytest.approx(final_flow)


def test_run_everything_v02_custom_input():
    custom_input_file = "../../test/input/yaml/CustomInput.yaml"
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

    run_parameters["nts"] = 12

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

    num_networks = len(results)
    id_list_j = 0
    values_list_j = 1
    file_base = "custom_test_standard_result"
    test_folder = Path("../../test/output/Pocono_TEST1/standard_test")

    results_standard = []
    for net_i in range(num_networks):
        ext_str_id = f".net{net_i:05d}_id.out" 
        ext_str_values = f".net{net_i:05d}_values.out" 
        with open(test_folder / (file_base + ext_str_id), "rb") as fid:
            ids = np.load(fid)
        with open(test_folder / (file_base + ext_str_values), "rb") as fid:
            values = np.load(fid)
        results_standard.append((ids,values))

        assert np.array_equal(results_standard[net_i][id_list_j], results[net_i][id_list_j])
        assert np.allclose(results_standard[net_i][values_list_j], results[net_i][values_list_j])


