import pathlib
from troute import HYFeaturesNetwork

hf_path = './test/hf/example.gpkg'

network = HYFeaturesNetwork.HYFeaturesNetwork(
    supernetwork_parameters={
        'geo_file_type': 'HYFeaturesNetwork',
        'geo_file_path': hf_path,
        'columns': {
            'key': 'id',
            'downstream': 'toid'
        }
    },
    waterbody_parameters={
        'break_network_at_waterbodies': False,
        'level_pool': {
            'level_pool_waterbody_parameter_file_path': hf_path,
            'reservoir_parameter_file': hf_path
        }
    },
    data_assimilation_parameters={},
    restart_parameters={
        'start_datetime': '2015-12-01 00:00:00'
    },
    compute_parameters={
        'parallel_compute_method': 'serial',
        'compute_kernel': 'V02-structured',
        'assume_short_ts': True,
        'cpu_pool': 6
    },
    forcing_parameters={
        'qlat_input_folder': './test/hf/qlat',
        'qts_subdivisions': 2,
        'dt': 2,
        'nts': 2
    },
    hybrid_parameters={},
    verbose=True
)

def _assemble():
    qlat_input_folder = pathlib.Path(network.forcing_parameters['qlat_input_folder'])
    all_files = sorted(qlat_input_folder.glob(network.forcing_parameters['qlat_file_pattern_filter']))
    assert(len(all_files) > 0)

    run_sets = network.build_forcing_sets()
    assert(len(run_sets) > 0)
    assert(len(run_sets[0]['qlat_files']) > 0)

    network.assemble_forcings(run_sets[0], )
    assert(network.qlateral.shape[0] > 0)

def test_read_csv():
    network.forcing_parameters['qlat_file_pattern_filter'] = "*NEXOUT.csv"
    _assemble()

def test_read_netcdf():
    network.forcing_parameters['qlat_file_pattern_filter'] = "*NEXOUT.nc"
    _assemble()
