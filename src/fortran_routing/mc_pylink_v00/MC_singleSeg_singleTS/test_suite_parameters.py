from itertools import zip_longest
import numpy as np

"""
see https://docs.google.com/spreadsheets/d/1kRdpoY2Ul0vZP2l9XXXT4tVA7QIPi8Co/edit#gid=1870807915
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
    qlat = (1, 500000)
    qup = (1, 500000)
    quc = (1, 500000)
    qdp = (1, 500000)
    depthp = (1, 100)
    dt = (5, 500)

    np.random.seed(randomseed)

    # rg = 15000

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
