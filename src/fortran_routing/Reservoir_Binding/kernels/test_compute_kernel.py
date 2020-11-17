"""
import pyximport
pyximport.install()
import sys
#Monkey patch pximport to not build modules not under test

finder = next(finder
              for finder in sys.meta_path
              if isinstance(finder, pyximport.PyxImporter))
find_module_original = finder.find_module

def find_module_new(self, fullname, package_path=None):
     if fullname.startswith('compute_kernel'):
         return find_module_original(fullname, package_path)

finder.find_module = find_module_new.__get__(finder, pyximport.PyxImporter)
"""
import pytest

from compute_kernel_lp import lp_kernel
from compute_kernel_hybrid import hybrid_kernel
from compute_kernel_rfc import rfc_kernel

@pytest.fixture()
def lp_reservoir():
    """
        an lp compute kernel to put under test
    """
    print ("lp test")

    water_elevation = 9.7373;
    lake_area = 15.0949;
    weir_elevation = 9.626;
    weir_coefficient = 0.4;
    weir_length = 10.0;
    dam_length = 10.0;
    orifice_elevation = 7.733;
    orifice_coefficient = 0.1;
    orifice_area = 1.0;
    max_depth = 9.96;
    lake_number = 16944276;

    k = lp_kernel(0, 1,
                water_elevation, lake_area,
                weir_elevation, weir_coefficient, weir_length,
                dam_length, orifice_elevation, orifice_coefficient,
                orifice_area, max_depth, lake_number)
    yield k


@pytest.fixture()
def hybrid_reservoir():
    """
        an lp compute kernel to put under test
    """
    print ("hybrid test")

    cwd_full = b"./reservoir_testing_files/";

    water_elevation = 1331.18005;
    lake_area = 209.632;
    weir_elevation = 1332.074;
    weir_coefficient = 0.4;
    weir_length = 10.0;
    dam_length = 10.0;
    orifice_elevation = 1314.473;
    orifice_coefficient = 0.1;
    orifice_area = 1.0;
    max_depth = 1335.180;
    initial_fractional_depth = 0.9
    lake_number = 402142;
    reservoir_type = 2;
    reservoir_parameter_file = b"./reservoir_testing_files/reservoir_index_short_range.nc";
    start_date = b"2010-10-01_07:00:00";
    usgs_timeslice_path = cwd_full;
    usace_timeslice_path = cwd_full;
    observation_lookback_hours = 48;
    observation_update_time_interval_seconds = 1000000000;

    k = hybrid_kernel(0, 1,
                water_elevation, lake_area,
                weir_elevation, weir_coefficient, weir_length,
                dam_length, orifice_elevation, orifice_coefficient,
                orifice_area, max_depth, initial_fractional_depth, 
                lake_number, reservoir_type, reservoir_parameter_file,
                start_date, usgs_timeslice_path, usace_timeslice_path, 
                observation_lookback_hours, observation_update_time_interval_seconds)
    yield k


@pytest.fixture()
def rfc_reservoir():
    """
        an lp compute kernel to put under test
    """
    print ("rfc test")

    cwd_full = b"./reservoir_testing_files/";

    water_elevation = 1331.18005;
    lake_area = 209.632;
    weir_elevation = 1332.074;
    weir_coefficient = 0.4;
    weir_length = 10.0;
    dam_length = 10.0;
    orifice_elevation = 1314.473;
    orifice_coefficient = 0.1;
    orifice_area = 1.0;
    max_depth = 1335.180;
    initial_fractional_depth = 0.9
    lake_number = 17609317;
    reservoir_type = 4;
    reservoir_parameter_file = b"./reservoir_testing_files/reservoir_index_short_range.nc";
    start_date = b"2019-08-18_09:00:00";
    time_series_path = cwd_full;
    forecast_lookback_hours = 24;


    k = rfc_kernel(0, 1,
                water_elevation, lake_area,
                weir_elevation, weir_coefficient, weir_length,
                dam_length, orifice_elevation, orifice_coefficient,
                orifice_area, max_depth, initial_fractional_depth, 
                lake_number, reservoir_type, reservoir_parameter_file,
                start_date, time_series_path, forecast_lookback_hours)
    yield k


def test_compute_kernel_lp(lp_reservoir):
    """
        test a simple construction of LP cython object
    """
    assert lp_reservoir is not None


def test_lp_run(lp_reservoir):
    """
        test running a LP reservoir
    """

    inflow_list = [91.27196, 91.7394, 92.15904, 92.1518, 91.84663, \
    91.38554, 90.86131, 90.32736, 89.81273, 89.3325, 88.89427, 88.5025, 88.16228, \
    87.41539, 86.80043, 86.03979, 85.3849, 85.33451, 86.84274, 91.6084, 101.81398, \
    118.85916, 143.99232, 177.7355, 219.2348, 267.22351, 319.90402, 374.54324, 428.86066, \
    480.92096, 529.23584, 572.77673, 610.93237, 643.4389, 670.28516, 691.67767, 707.96088, \
    719.57312, 726.96997, 730.63269, 731.03186, 728.61438, 723.79578, 716.9549, 708.43268, \
    698.53247, 687.52112, 675.63123, 663.06421, 649.99976, 636.57898, 622.92926, 609.1745, \
    595.40369, 581.68799, 568.08588, 554.64484, 541.4032, 528.39185, 515.63513, 503.14838, \
    490.95123, 479.05109, 467.45493, 456.16663, 445.18753, 434.51706, 424.15311,414.0921, \
    404.32956, 394.86014, 385.67789, 376.77621, 368.14966, 359.78958, 351.68875, 343.83972, \
    336.23505, 328.86719, 321.7287, 314.81219, 308.11047, 301.61646, 295.32312, 289.22369, \
    283.31207, 277.5813, 272.02521, 266.63776, 261.41315, 256.34564, 251.42978, 246.66023, \
    242.03192, 237.53989, 233.17944, 228.94595, 224.83511, 220.84265, 216.96449, 213.19672, \
    209.53554, 205.97734, 202.51857, 199.1559, 195.88605, 192.70595, 189.61255]

    routing_period = 300.0

    for inflow in inflow_list:
        out = lp_reservoir.run(inflow, 0.0, routing_period)
        print(out)

    expected_final_outflow = 17.0437641

    assert expected_final_outflow == pytest.approx(out)


def test_compute_kernel_hybrid(hybrid_reservoir):
    """
        test a simple construction of hybrid cython object
    """
    assert hybrid_reservoir is not None


def test_compute_hybrid_run(hybrid_reservoir):
    """
        test running a hybrid reservoir
    """

    inflow_list = [189.22899, 189.27005, 189.31049, 189.35042, 189.38965, 189.42819, 189.46588, 189.50273, \
    189.53859, 189.57346, 189.60719, 189.63979, 189.6711, 189.7011, 189.72968, \
    189.75679, 189.7823, 189.80617, 189.82822, 189.84842, 189.86653, 189.88255, 189.89622, \
    189.90752, 189.91612, 189.922, 189.92482, 189.92447, 189.92067, 189.91319, 189.90175, \
    189.88611, 189.86592, 189.84088, 189.81064, 189.77487, 189.73317, 189.6852, 189.63051, \
    189.56873, 189.49939, 189.42207, 189.33635, 189.24176, 189.13782, 189.02408, \
    188.90009, 188.76535, 188.61945, 188.46188, 188.29224, 188.11006, 187.91493, 187.70644, \
    187.48419, 187.24779, 186.9969, 186.73119, 186.45035, 186.15407, 185.84213, 185.51424, \
    185.17023, 184.80989, 184.43312, 184.03975, 183.62973, 183.20296, 182.75943, 182.29909, \
    181.82205, 181.32828, 180.81792, 80.29099, 179.74774, 179.1882, 178.61267, 178.02129, \
    177.41437, 176.79207, 176.15475, 175.50269, 174.83627, 174.15576, 173.46162, \
    172.75417, 172.03389, 171.3011, 170.55634, 169.79997, 169.03255, 168.25441, 167.46616, \
    166.66815, 165.86099, 165.04509, 164.22101, 163.38913, 162.55011, 161.70428, 160.85229, \
    159.99452, 159.13156, 158.26382, 157.39188, 156.51611, 155.63715, 154.75531, 153.8712, 152.98517, \
    152.09779, 151.2094, 150.32057, 149.43166, 148.54315, 147.6554, 146.76892, 145.88405, 145.00128, 144.12091]

    routing_period = 300.0

    for inflow in inflow_list:
        out = hybrid_reservoir.run(inflow, 0.0, routing_period)
        print(out)

    expected_final_outflow = 13.73367

    assert expected_final_outflow == pytest.approx(out)


def test_compute_kernel_rfc(rfc_reservoir):
    """
        test a simple construction of RFC cython object
    """
    assert rfc_reservoir is not None


def test_compute_rfc_run(rfc_reservoir):
    """
        test running a RFC reservoir
    """

    inflow_list = [189.22899, 189.27005, 189.31049, 189.35042, 189.38965, 189.42819, 189.46588, 189.50273, \
    189.53859, 189.57346, 189.60719, 189.63979, 189.6711, 189.7011, 189.72968, \
    189.75679, 189.7823, 189.80617, 189.82822, 189.84842, 189.86653, 189.88255, 189.89622, \
    189.90752, 189.91612, 189.922, 189.92482, 189.92447, 189.92067, 189.91319, 189.90175, \
    189.88611, 189.86592, 189.84088, 189.81064, 189.77487, 189.73317, 189.6852, 189.63051, \
    189.56873, 189.49939, 189.42207, 189.33635, 189.24176, 189.13782, 189.02408, \
    188.90009, 188.76535, 188.61945, 188.46188, 188.29224, 188.11006, 187.91493, 187.70644, \
    187.48419, 187.24779, 186.9969, 186.73119, 186.45035, 186.15407, 185.84213, 185.51424, \
    185.17023, 184.80989, 184.43312, 184.03975, 183.62973, 183.20296, 182.75943, 182.29909, \
    181.82205, 181.32828, 180.81792, 80.29099, 179.74774, 179.1882, 178.61267, 178.02129, \
    177.41437, 176.79207]

    routing_period = 3600.0

    for inflow in inflow_list:
        out = rfc_reservoir.run(inflow, 0.0, routing_period)
        print(out)

    expected_final_outflow = 3.6

    assert expected_final_outflow == pytest.approx(out)

