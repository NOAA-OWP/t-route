import traceback

debuglevel = -1
COMPILE = True
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r"f2py3")
        fortran_compile_call.append(r"-c")
        fortran_compile_call.append(r"MCsingleSegStime_f2py_NOLOOP.f90")
        fortran_compile_call.append(r"-m")
        fortran_compile_call.append(r"mc_sseg_stime")
        if debuglevel <= -2:
            subprocess.run(fortran_compile_call)
        else:
            subprocess.run(
                fortran_compile_call,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        import mc_sseg_stime as mc
    except Exception as e:
        print(e)
        if debuglevel <= -1:
            traceback.print_exc()
else:
    import mc_sseg_stime as mc

# # Method 1
# Python: time loop; segment loop; constant channel variables are passed to Fortran
# Fortran: Take constant variable values and then run MC for a single segment


def compute_mc_up2down_ReachbySegment():
    """HANDLE LOOPING, 
        Then call single segment routine for each segment"""
    pass


def singlesegment(
    dt,  # dt
    qup=None,  # qup
    quc=None,  # quc
    qdp=None,  # qdp
    qlat=None,  # ql
    dx=None,  # dx
    bw=None,  # bw
    tw=None,  # tw
    twcc=None,  # twcc
    n_manning=None,  #
    n_manning_cc=None,  # ncc
    cs=None,  # cs
    s0=None,  # s0
    velp=None,  # velocity at previous time step
    depthp=None,  # depth at previous time step
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
    # return qdc, vel, depth


def main():

    precision = "double"  # the fortran real*8 declarations must be changed to
    # precision = 'single # real*4 to use precision single and obtain the expected output below

    dt = 60.0  # Time step
    dx = 1800.0  # segment length
    bw = 112.0  # Trapezoidal bottom width
    tw = 448.0  # Channel top width (at bankfull)
    twcc = 623.5999755859375  # Flood plain width
    n_manning = 0.02800000086426735  # manning roughness of channel
    n_manning_cc = 0.03136000037193298  # manning roughness of floodplain
    cs = 1.399999976158142  # channel trapezoidal sideslope
    s0 = 0.0017999999690800905  # downstream segment bed slope
    qlat = 40.0  # Lateral inflow in this time step

    if precision is "single":
        """
        single precision results from standard procedure
        k,   i,   q,    vel,    depth
        0 0 0.0891760066151619 0.06941037625074387 0.010050795041024685
        0 1 0.06863606721162796 0.07104580104351044 0.009998836554586887
        0 2 0.04598825052380562 0.06900732219219208 0.009981763549149036
        0 3 0.21487340331077576 0.07048019021749496 0.010033470578491688
        1 0 0.3173820674419403 0.15248453617095947 0.03275299444794655
        1 1 0.2014089822769165 0.09322302788496017 0.015030190348625183
        1 2 0.13257014751434326 0.061184611171483994 0.008333344012498856
        1 3 0.7570106983184814 0.12373604625463486 0.02334451675415039
        """
        qup = 0.04598825052380562  # Flow from the upstream neighbor in the previous timestep
        quc = 0.04598825052380562  # Flow from the upstream neighbor in the current timestep
        qdp = (
            0.21487340331077576  # Flow at the current segment in the previous timestep
        )
        depthp = 0.010033470578491688  # Depth at the current segment in the previous timestep')
        velp = 0.07048019021749496  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!

    elif precision is "double":
        """
        double precision results from standard procedure
        k,   i,   q,    vel,    depth
        0 0 0.08917599589412113 0.06941036720349936 0.010050795413704227
        0 1 0.06863606840763538 0.07104580614454808 0.009998836839100517
        0 2 0.04598825231793033 0.06900732835665264 0.009981763746485387
        0 3 0.21487341391047937 0.07048019432213874 0.01003347079559232
        1 0 0.3173820852675605 0.15248449837070724 0.03275299172786703
        1 1 0.20140897858530027 0.09322303229159627 0.01503019048616059
        1 2 0.13257014557875424 0.061184620519535104 0.008333345170456535
        1 3 0.7570106991484074 0.12373605347605196 0.02334451636521022
        """
        qup = 0.04598825231793033  # Flow from the upstream neighbor in the previous timestep
        quc = 0.04598825231793033  # Flow from the upstream neighbor in the current timestep
        qdp = (
            0.21487341391047937  # Flow at the current segment in the previous timestep
        )
        depthp = 0.01003347079559232  # Depth at the current segment in the previous timestep')
        velp = 0.07048019432213874  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!

    qdc_expected = 0.7570106983184814
    velc_expected = 0.12373604625463486
    depthc_expected = 0.02334451675415039

    qdc_expected_dbl = 0.7570106991484074
    velc_expected_dbl = 0.12373605347605196
    depthc_expected_dbl = 0.02334451636521022

    # run M-C model
    qdc, velc, depthc = singlesegment(
        dt=dt,
        qup=qup,
        quc=quc,
        qdp=qdp,
        qlat=qlat,
        dx=dx,
        bw=bw,
        tw=tw,
        twcc=twcc,
        n_manning=n_manning,
        n_manning_cc=n_manning_cc,
        cs=cs,
        s0=s0,
        velp=velp,
        depthp=depthp,
    )
    print(
        "real double precision computed q: {} vel: {} depth: {}".format(
            qdc, velc, depthc
        )
    )
    print(
        "real double precision expected q: {} vel: {} depth: {}".format(
            qdc_expected_dbl, velc_expected_dbl, depthc_expected_dbl
        )
    )
    print(
        "real single precision expected q: {} vel: {} depth: {}".format(
            qdc_expected, velc_expected, depthc_expected
        )
    )


if __name__ == "__main__":
    main()
