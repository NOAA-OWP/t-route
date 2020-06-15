import sys

sys.path.append(r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS")
# import mc_sseg_stime_NOLOOP as mc
from mc_sseg_stime import muskingcunge_module as mc


def singlesegment(
    dt,  # dt
    qup=None,  # qup can be anything
    quc=None,  # quc will not be more than a 10 percent diff than qup
    qlat=None,  # ql can be anything - key
    qdp=None,  # qdp will not be more than 20 percent diff than qup+qlat
    dx=None,  # dx fully variable
    bw=None,  # bw correlated to tw, tw always > bw
    tw=None,  # tw correlated to bw, bw always < tw
    twcc=None,  # twcc always > than tw, tw of broader floodplain
    n_manning=None,  # usually smaller than n_manning_cc
    n_manning_cc=None,  # ncc usually greater than n_manning
    cs=None,  # cs correlated to bw and tw
    s0=None,  # s0 variable
    velp=None,  # velocity at previous time step not rel
    depthp=None,  # depth at previous time step starting point for iteration depthp = approx(y_direct(bw,n_manning,s0,avg(qup,qdp)))
):

    # call Fortran routine
    return mc.muskingcungenwm(
        dt,
        qup,
        quc,
        qdp,
        qlat,
        dx,
        bw,
        tw,
        twcc,
        n_manning,
        n_manning_cc,
        cs,
        s0,
        velp,
        depthp,
    )
