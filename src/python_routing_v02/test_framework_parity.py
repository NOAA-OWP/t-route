from nwm_routing.__main__ import _run_everything_v02  # , _input_handler_v02
# from nwm_routing.__main__ import _run_everything_v03  # , _input_handler_v03
# from troute.routing.compute import compute_nhd_routing_v02
from troute.nhd_io import read_custom_input
from pathlib import Path
import numpy as np
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


