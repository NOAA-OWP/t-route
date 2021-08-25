import pytest
import mc_sseg_stime_NOLOOP_demo as demo
from test_suite_parameters import generate_conus_MC_parameters


"""
Statistics from the CONUS NWM 2.1 dataset
            mean     std_dev      min        25%         50%          75%           max
dx    1,947.776    1,965.625     1.000   558.000   1,549.000    2,631.000    95,714.000 
bw        5.290        9.407     0.135     1.904       2.769        4.889       230.035 
tw        8.817       15.678     0.225     3.173       4.615        8.148       383.392 
twcc     26.450       47.033     0.674     9.518       13.845      24.443     1,150.175 
n         0.058        0.003     0.040     0.060       0.060        0.060         0.060 
ncc       0.117        0.006     0.080     0.120       0.120        0.120         0.120 
cs        0.5857       0.1945    0.0846    0.4625      0.5944       0.7012        2.254 
s0        0.02150      0.04585   0.00001   0.00100     0.00600      0.01900       4.600 
qlat                              1       500000
qup                               1       500000
quc                               1       500000
qdp                               1       500000
depthp                            1       100
dt==300

*Legend*
dt = Time step
dx = segment length
bw = Trapezoidal bottom width
tw = Channel top width (at bankfull)
twcc = Flood plain width
n_manning = manning roughness of channel
n_manning_cc = manning roughness of floodplain
cs = channel trapezoidal sideslope
s0 = downstream segment bed slope
qlat = Lateral inflow in this time step
qup = inflow into segment from upstream in previous timestep
quc = inflow into segment from upstream in current timestep
qdp = outflow from segment in previous timestep
depthp = depth in segment in previous timestep
qdc = outflow from segment in current timestep (MC solution yields this value)

*Notes:*
ncc == 2 * n

twcc == 3 * tw

tw >= bw

prior tests have varied dt, even though it is always held at 300 for NWM simulations to date.
"""

input_tup = (
    "single",  # calculation precision (NOT USED)
    300.0,  # Time step
    1800.0,  # segment length
    112.0,  # Trapezoidal bottom width
    448.0,  # Channel top width (at bankfull)
    623.5999755859375,  # Flood plain width
    0.02800000086426735,  # manning roughness of channel
    0.03136000037193298,  # manning roughness of floodplain
    1.399999976158142,  # channel trapezoidal sideslope
    0.0017999999690800905,  # downstream segment bed slope
    40.0,  # Lateral inflow in this time step
    4509,
    5098,
    5017,
    30,
)

# qdc1, qdc2, velc1, velc2, depthc1, depthc2 = demo.compare_methods(*input_tup)

input_tup_random = generate_conus_MC_parameters(3, 10)
out = demo.compare_methods("single", *(next(input_tup_random)))
# print(out)
qdc1, qdc2, velc1, velc2, depthc1, depthc2 = out


def test_MC_kernel_q():
# TODO: Take advantage of ranges above to build a state-space 
# exploration of potential inputs and confirm it parity across all inputs

    assert qdc1 == qdc2


def test_MC_kernel_vel():
# TODO: Take advantage of ranges above to build a state-space 
# exploration of potential inputs and confirm it parity across all inputs

    assert velc1 == velc2
 

def test_MC_kernel_depth():
# TODO: Take advantage of ranges above to build a state-space 
# exploration of potential inputs and confirm it parity across all inputs

   assert depthc1 == depthc2


