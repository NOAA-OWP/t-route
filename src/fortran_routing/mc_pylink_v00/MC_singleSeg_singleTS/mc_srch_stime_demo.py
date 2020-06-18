from __future__ import absolute_import
from __future__ import print_function
import os
import traceback
import subprocess
import numpy as np

debuglevel = 0
COMPILE = True
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r"f2py3")
        fortran_compile_call.append(r"-c")
        fortran_compile_call.append(r"varPrecision.f90")
        fortran_compile_call.append(r"varSingleSegStime_f2py_2clean.f90")
        fortran_compile_call.append(r"MCsingleSegStime_f2py_2clean.f90")
        fortran_compile_call.append(r"-m")
        fortran_compile_call.append(r"mc_srch_stime")
        # fortran_compile_call.append(r"--opt='-fdefault-real-8'")
        if debuglevel <= -2:
            subprocess.run(fortran_compile_call)
        else:
            subprocess.run(
                fortran_compile_call,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        import mc_srch_stime as mc
    except Exception as e:
        print(e)
        if debuglevel <= -1:
            traceback.print_exc()
else:
    import mc_srch_stime as mc


# # Method 1
# Python: time loop; segment loop; constant channel variables are passed to Fortran
# Fortran: Take constant variable values and then run MC for a single segment

# # test input dataset

# In[6]:


ntim = 2
# the number of time steps necessary for variables passed to mc module to compute successfully
nlinks = 2
# the number of links needed to define varialbe qd. ** nlinks is not used in fortran source code.
mc.var.uslinkid = 1
mc.var.linkid = 2

ncomp0 = 2
mc.var.ncomp0 = (
    ncomp0  # the number of segments of a reach upstream of the current reach
)
ncomp = 4
mc.var.ncomp = ncomp  # the number of segments of the current reach
mxseg = max(ncomp0, ncomp)


# MC model key arrays
# mc.var.qd=np.zeros((ntim,mxseg,nlinks))  #will store MC output qdc (flow downstream current timestep)
mc.var.qd = np.zeros(
    (ntim, ncomp)
)  # will store MC output qdc (flow downstream current timestep)
mc.var.vela = np.zeros((ntim, ncomp))
mc.var.deptha = np.zeros((ntim, ncomp))


# dt
dtp = np.zeros((ncomp))
dtp[0] = 60.0
dtp[1] = 40.0
dtp[2] = 20.0
dtp[3] = 60.0
# dx
dxp = np.zeros((ncomp))
dxp[0] = 1200.0
dxp[1] = 2500.0
dxp[2] = 2800.0
dxp[3] = 1800.0
# bw
bwp = np.zeros((ncomp))
bwp[0] = 50.0
bwp[1] = 100.0
bwp[2] = 148.0
bwp[3] = 112.0
# tw
twp = np.zeros((ncomp))
twp[0] = 4.0 * bwp[0]
twp[1] = 4.2 * bwp[1]
twp[2] = 4.8 * bwp[2]
twp[3] = 4.0 * bwp[3]
# twcc
twccp = np.zeros((ncomp))
twccp[0] = twp[0] + 103.2
twccp[1] = twp[1] + 204.4
twccp[2] = twp[2] + 275.3
twccp[3] = twp[3] + 175.6
# n
npy = np.zeros((ncomp))
npy[0] = 0.03
npy[1] = 0.032
npy[2] = 0.038
npy[3] = 0.028
# ncc
nccp = np.zeros((ncomp))
nccp[0] = 1.1 * npy[0]
nccp[1] = 1.2 * npy[1]
nccp[2] = 1.24 * npy[2]
nccp[3] = 1.12 * npy[3]
# cs
csp = np.zeros((ncomp))
csp[0] = 1.2
csp[1] = 1.8
csp[2] = 2.6
csp[3] = 1.4
# so
sop = np.zeros((ncomp))
sop[0] = 0.002
sop[1] = 0.0024
sop[2] = 0.0032
sop[3] = 0.0018
# lateral flow
qlat = np.zeros((ncomp))


# In[3]:


#%%pixie_debugger

# *** WARNING: All Python variables corresponding to Fortran module must be in lowercase.
# fortran source code that is used for f2py transformation includes do k=1, ntim loop so
# that multiple timestep computation is activated inside of the f2py module.

# run M-C model
nts = 2  # the number of timestep in simulation
# variable storing all outputs in time
wqd = np.zeros((nts, mxseg))
wvela = np.zeros((nts, mxseg))
wdeptha = np.zeros((nts, mxseg))
# TODO: Add a few more diagnostic outputs -- this array could be populated and printed, for instance.
qla = np.zeros((nts, mxseg))

bply = nts / 2.0
# uslinkflag=1/0 when upstream reach exist/not exist:
mc.var.uslinkflag = 1

for k in range(0, nts):
    # upstream reach one timestep behind flow. Note that for qd[k,i,j]
    # k=0/1: previous/current timestep; i: node ID; j=0/1: upstream/current reach
    aply = -0.2
    # mc.var.qd[0,ncomp0-1,0]= aply*(k - bply)**2.0 - aply*bply**2.0
    mc.var.qus_prev = aply * (k - bply) ** 2.0 - aply * bply ** 2.0

    for i in range(0, ncomp):
        # lateral flow for current reach.
        aply = -0.8
        qlat[i] = aply * (k + 1 - bply) ** 2.0 - aply * bply ** 2.0 + (i + 1) * 10.0
        mc.var.ql = qlat[i]
        # channel data
        mc.var.dt = dtp[i]
        mc.var.dx = dxp[i]
        mc.var.bw = bwp[i]
        mc.var.tw = twp[i]
        mc.var.twcc = twccp[i]
        mc.var.n = npy[i]
        mc.var.ncc = nccp[i]
        mc.var.cs = csp[i]
        mc.var.so = sop[i]
        # current node
        mc.var.iseg = i + 1
        # initial values
        if k == 0:
            # initialize upstream reach
            # mc.var.qd[0,ncomp0-1,0]= 0.0
            mc.var.qus_prev = 0.0
            # initialize current reach
            for ii in range(0, ncomp):
                # mc.var.qd[0,i,1]= 0.0
                mc.var.qd[0, ii] = 0.0
                mc.var.vela[0, ii] = 0.0
                mc.var.deptha[0, ii] = -999.0

        # call Fortran routine
        mc.muskingcunge_module.main()

        # print channel data
        debuglevel = -2
        if debuglevel <= -2:
            str = ""
            str = str + f" dt: {mc.var.dt}"
            str = str + f" dx: {mc.var.dx}"
            str = str + f" bw: {mc.var.bw}"
            str = str + f" tw: {mc.var.tw}"
            str = str + f" twcc: {mc.var.twcc}"
            str = str + f" n: {mc.var.n}"
            str = str + f" ncc: {mc.var.ncc}"
            str = str + f" cs: {mc.var.cs}"
            str = str + f" so: {mc.var.so}"
            str = str + f" ql: {mc.var.ql}"
            str = str + f" qup: {mc.var.qup}"
            str = str + f" quc: {mc.var.quc}"
            str = str + f" qdp: {mc.var.qdp}"

            str = str + f" velp_check: {mc.var.velp_chk}"
            str = str + f" depthp_check: {mc.var.depthp_chk}"
            str = str + f" qdc: {mc.var.qdc}"
            str = str + f" vel: {mc.var.vel}"
            str = str + f" depth: {mc.var.depth}"

            print(str)

        # store computed values at time k for computation at k+1
        # qd[k,i,j]: k=0/1: previous/current timestep; i: node ID; j=0/1: upstream/current reach
        # mc.var.qd[0,i,1]= mc.var.qd[1,i,1]
        # output keeping
        wqd[k, i] = mc.var.qd[1, i]
        wvela[k, i] = mc.var.vela[1, i]
        wdeptha[k, i] = mc.var.deptha[1, i]

    for i in range(0, ncomp):
        mc.var.qd[0, i] = wqd[k, i]
        mc.var.vela[0, i] = wvela[k, i]
        mc.var.deptha[0, i] = wdeptha[k, i]


# test output
# =1
print(r"k,   i,   q,    vel,    depth")
for k in range(0, nts):
    for i in range(0, ncomp):
        print(k, i, wqd[k, i], wvela[k, i], wdeptha[k, i])
