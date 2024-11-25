from itertools import zip_longest
import numpy as np

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

These values are bounded by the estimated maximum discharge (70k cms) at
Arkansas City during the 1927 flood on the Mississippi.
See: https://pubs.usgs.gov/circ/2004/circ1254/pdf/circ1254.pdf
The deepest river in the US is the Hudson, at around 200 ft. The model
may often seed greater depths, because of the geometric formulation of
the current NWM, so the test range is given a buffer.
        min      max
qlat                              1       70000
qup                               1       70000
quc                               1       70000
qdp                               1       70000
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

def generate_conus_MC_parameters(rg = 1, randomseed = None):
    dx = ( 1947.776, 1965.625, 1.000, 558.000, 1549.000, 2631.000, 95714.000,)
    bw = ( 5.290, 9.407, 0.135, 1.904, 2.769, 4.889, 230.035,)
    tw = ( 8.817, 15.678, 0.225, 3.173, 4.615, 8.148, 383.392,)
    twcc = ( 26.450, 47.033, 0.674, 9.518, 13.845, 24.443, 1150.175,)
    n = ( 0.058, 0.003, 0.040, 0.060, 0.060, 0.060, 0.060,)
    ncc = ( 0.117, 0.006, 0.080, 0.120, 0.120, 0.120, 0.120,)
    cs = ( 0.5857, 0.1945, 0.0846, 0.4625, 0.5944, 0.7012, 2.254,)
    s0 = ( 0.02150, 0.04585, 0.00001, 0.00100, 0.00600, 0.01900, 4.600,)
    qlat = (1, 70000)
    qup = (1, 70000)
    quc = (1, 70000)
    qdp = (1, 70000)
    depthp = (1, 25)
    dt = (5, 500)

    np.random.seed(randomseed)

    test_suite_parameter_set = zip_longest(
        np.random.uniform(dx[2], dx[6], rg),
        np.random.uniform(bw[2], dx[6], rg),
        np.random.uniform(tw[2], tw[6], rg),
        np.random.uniform(twcc[2], twcc[6], rg),
        np.random.uniform(n[2], n[6], rg),
        np.random.uniform(ncc[2], ncc[6], rg),
        np.random.uniform(cs[2], cs[6], rg),
        np.random.uniform(s0[2], s0[6], rg),
        np.random.uniform(qlat[0], qlat[1], rg),
        np.random.uniform(qup[0], qup[1], rg),
        np.random.uniform(quc[0], quc[1], rg),
        np.random.uniform(qdp[0], qdp[1], rg),
        np.random.uniform(depthp[0], depthp[1], rg),
        np.random.uniform(dt[0], dt[1], rg),
    )

    # print("\n".join([f"{t}" for t in test_suite_parameter_set]))

    return test_suite_parameter_set
