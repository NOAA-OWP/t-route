import traceback

debuglevel = 0
COMPILE = True
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r"f2py3")
        fortran_compile_call.append(r"-c")
        # fortran_compile_call.append(r"varPrecision.f90")
        fortran_compile_call.append(r"module_levelpool.f90")
        fortran_compile_call.append(r"-m")
        fortran_compile_call.append(r"levelpools")
        # fortran_compile_call.append(r"--opt='-fdefault-real-8'")

        if debuglevel <= -1:
            print(fortran_compile_call)
        if debuglevel <= -2:
            subprocess.run(fortran_compile_call)
        else:
            subprocess.run(
                fortran_compile_call,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        from levelpools import *

        # import mc_sseg_stime as muskingcunge_module
    except Exception as e:
        print(e)
        if debuglevel <= -1:
            traceback.print_exc()
else:
    from levelpools import *

    # import mc_sseg_stime as muskingcunge_module

# # Method 1
# Python: time loop; segment loop; constant channel variables are passed to Fortran
# Fortran: Take constant variable values and then run MC for a single segment

#comment this section out when attaching new data, fortran code needs values to run correctly
# ln=1,
# qi0=1,
# qi1=1,
# ql=1,
# dt=1,
# h=1,
# ar=1,
# we=1,
# maxh=1,
# wc=1,
# wl=1,
# dl=1,
# oe=1,
# oc=1,
# oa=1,


def reservoirs_calc(
    ln=None,
    qi0=None,
    qi1=None,
    ql=None,
    dt=None,
    h=None,
    ar=None,
    we=None,
    maxh=None,
    wc=None,
    wl=None,
    dl=None,
    oe=None,
    oc=None,
    oa=None,
):

    # call Fortran routine
    return module_levelpool.levelpool_physics(
    ln,
    qi0,
    qi1,
    ql,
    dt,
    h,
    ar,
    we,
    maxh,
    wc,
    wl,
    dl,
    oe,
    oc,
    oa,
    )
    # return qdc, vel, depth


# reservoirs_calc()
def main():
    reservoirs_calc(
    ln,
    qi0,
    qi1,
    ql,
    dt,
    h,
    ar,
    we,
    maxh,
    wc,
    wl,
    dl,
    oe,
    oc,
    oa,
)



if __name__ == "__main__":
    main()
