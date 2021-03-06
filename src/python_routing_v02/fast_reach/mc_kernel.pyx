# cython: language_level=3, cdivision=True

cdef extern from "<math.h>" nogil:
    double fabs(double x)
    float fabsf(float)
    double fmin(double x, double y)
    float fminf(float, float)
    double fmax(double x, double y)
    float fmaxf(float, float)
    double sqrt(double x)
    float sqrtf(float)
    double pow(double x, double y)
    float powf(float, float)

cdef struct QC:
    float Qj
    float Qj_0
    float C1
    float C2
    float C3
    float C4


# INPUTS:
# Arrays:
# Q: np.array([qup, quc, qdp, 1])
# C: np.array([C1, C2, C3, C4])
# dx, bw, tw, twcc, ncc, cs, so, n, qlat: 1D arrays
# vela, deptha: 2D arrays
# qd: 3D arrays

# Scalars:
# dt, qup, quc, qdp, qdc, ql, z, vel, depth: float
# bfd, WPC, AREAC: float
# ntim: integer
# ncomp, linkID: integer

# dt = 60.0 # Time step
# dx = 1800.0 # segment length
# bw = 112.0 # Trapezoidal bottom width
# tw = 448.0 # Channel top width (at bankfull)
# twcc = 623.5999755859375 # Flood plain width
# n = 0.02800000086426735 # manning roughness of channel
# ncc = 0.03136000037193298 # manning roughness of floodplain
# cs = 1.399999976158142 # channel trapezoidal sideslope
# s0 = 0.0017999999690800905 # downstream segment bed slope
# ql = 40.0 # Lateral inflow in this time step
# qup = 0.04598825052380562 # Flow from the upstream neighbor in the previous timestep
# quc = 0.04598825052380562 # Flow from the upstream neighbor in the current timestep
# qdp = 0.21487340331077576 # Flow at the current segment in the previous timestep
# depthp = 0.010033470578491688 # Depth at the current segment in the previous timestep')
# Expected (0.7570107129107292, 0.12373605608067537, 0.023344515869393945)

#args = (60.0, 0.04598825052380562, 0.04598825052380562, 0.21487340331077576,
 #       40.0, 1800.0, 112.0, 448.0, 623.5999755859375, 0.02800000086426735,
 #       0.03136000037193298, 1.399999976158142, 0.0017999999690800905, 0.010033470578491688)
#
#single_precision = tuple(map(np.float32, args))
#double_precision = tuple(map(np.float64, args))

# single_seg.muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc, n, ncc, cs, s0, depthp)

cdef void cython_muskingcunge(
    float dt,
		float qup,
		float quc,
		float qdp,
		float ql,
		float dx,
		float bw,
		float tw,
		float twcc,
		float n,
		float ncc,
		float cs,
		float s0,
		float velp,
		float depthp,
    QVD *rv,
) nogil:
    cdef int maxiter = 100
    cdef float mindepth = 0.01
    cdef float aerror = 0.01
    cdef float rerror = 1.0
    cdef int it, tries = 0

    cdef float z, h, h_0
    cdef float qdc, velc, depthc
    cdef float C_pdot_Q
    cdef float R, twl
    cdef float bfd
    cdef QC qc_struct
    # initialize vars
    qc_struct.Qj_0 = 0.0
    qc_struct.Qj = 0.0
    qc_struct.C1 = 0.0
    qc_struct.C2 = 0.0
    qc_struct.C3 = 0.0
    qc_struct.C4 = 0.0

    cdef QC *qc = &qc_struct

    cdef float C1, C2, C3, C4
    cdef float Qj_0, Qj

    if cs == 0:  # susceptible to float comparison error?
        z = 1
    else:
        z = 1/cs

    if bw > tw:
        bfd = bw * (1/0.00001)
    elif bw == tw:
        bfd = bw/(2*z)
    else:
        bfd = (tw - bw)/(2*z)

    depthc = fmaxf(depthp, 0)
    h = (depthc * 1.33) + mindepth # 4/3 ?
    h_0 = depthc * 0.67  # 2/3 ?

    if ql > 0 or quc > 0 or qup > 0 or qdp > 0:

        # These should not be carried between iterations --- see https://github.com/NCAR/wrf_hydro_nwm_public/pull/510
        # WPC = 0
        # AREAC = 0

        it = 0

        while rerror > 0.01 and aerror >= mindepth and it <= maxiter:
            secant_h0(z, bw, bfd, twcc, s0, n, ncc, dt, 
                      dx, qup, quc, qdp, ql, h_0,
                      qc)
            secant_h(z, bw, bfd, twcc, s0, n, ncc, dt, 
                      dx, qup, quc, qdp, ql, h, 
                      qc)

            Qj_0 = qc.Qj_0
            Qj = qc.Qj

            C1 = qc_struct.C1
            C2 = qc_struct.C2
            C3 = qc_struct.C3
            C4 = qc_struct.C4

            if Qj_0 - Qj != 0:
                h_1 = h - (Qj * (h_0 - h)) / (Qj_0 - Qj)

                if h_1 < 0:
                    h_1 = h
            else:
                h_1 = h

            if h > 0:
                rerror = fabsf((h_1 - h) * (1/h))
                aerror = fabsf(h_1 - h)
            else:
                rerror = 0
                aerror = 0.9

            h_0 = max(0, h)
            h = max(0, h_1)
            it += 1

            if h < mindepth:
                if it >= maxiter:
                    tries += 1
                    if tries <= 4:
                        h = h * 1.33
                        h_0 = h_0 * 0.67
                        maxiter += 25
                        Qj_0 = 0
                        WPC = 0
                        AREAC = 0
                        _iter = 0
                        continue
                    else:
                        break
        C_pdot_Q = (C1 * qup) + (C2 * quc) + (C3 * qdp)
        qdc = C_pdot_Q + C4
        if qdc < 0:
            if C4 < 0 and fabsf(C4) > C_pdot_Q:
                qdc = 0
            else:
                qdc = max((C1 * qup) * (C2 * quc) + C4, (C1 * qup) + (C3 * qdp) + C4)

        twl = bw + (2 * z * h)
        R = (0.5 * h * (bw + twl)) / (bw + 2.0 * sqrtf(powf(((twl - bw) * 0.5), 2.0) + (h * h)))
        velc = (1 / n) * (powf(R, 2.0/3.0)) * sqrtf(s0)
        depthc = h
    else:
        qdc = 0
        depthc = 0
        velc = 0

    rv.qdc = qdc
    rv.velc = velc
    rv.depthc = depthc


cdef inline void secant_h0(
    float z,
		float bw,
		float bfd,
		float twcc,
		float s0,
		float n,
		float ncc,
		float dt,
		float dx,
		float qup,
		float quc,
		float qdp,
		float ql,
		float h_0,
		QC* qc,
) nogil:
    """
    Returns Qj_0, C
    """
    cdef float twl, AREA, AREAC, WP, WPC, R, Ck, X, Km, D

    # top surface water width of the channel inflow
    twl = bw + 2*z*h_0

    cdef float sqr_s0 = sqrtf(s0)
    cdef float sqr_1z2 = sqrtf(1.0 + (z * z))
    cdef float areasum

    if h_0 > bfd:
        # hydraulic radius, R
        AREA = (bw+bfd*z) * bfd
        AREAC = twcc * (h_0 - bfd)
        WP = bw + 2.0 * bfd * sqrtf(1 + (z * z))
        WPC = twcc + (2 * (h_0 - bfd))
        R = (AREA + AREAC)/(WP + WPC)

        # kinematic celerity, c
        areasum = (1/(AREA + AREAC))
        Ck = fmaxf(0.0, ((sqr_s0 / n) * ((5./3.) * powf(R,(2./3.)) -
                         ((2./3.) * powf(R,(5. / 3.)) * (2.0 * sqr_1z2 / (bw + 2.0 * bfd * z)))) * AREA
                       + ((sqr_s0 / ncc) * (5. / 3.) * powf((h_0 - bfd),(2. / 3.))) * AREAC) * areasum)

        # MC parameter, X
        X = fminf(0.5, fmaxf(0.0, 0.5 * (1 - (0.5 * (qc.Qj_0 / (twcc * s0 * Ck * dx))))))
    else:
        # hydraulic radius, R
        AREA = (bw + h_0 * z) * h_0
        AREAC = 0
        WP = bw + 2 * h_0 * sqr_1z2
        WPC = 0
        if WP > 0:
            R = AREA / WP
        else:
            R = 0

        # kinematic celerity, c
        if h_0 > 0:
            Ck = fmaxf(0.0, (sqr_s0 / n) * ((5. / 3.) * powf(R,(2. / 3.)) -
                                        ((2. / 3.) * powf(R,(5. / 3.)) * (
                                                    2.0 * sqr_1z2 / (bw + 2.0 * h_0 * z)))))
            X = fminf(0.5, fmaxf(0.0, 0.5 * (1 - (0.5 * (qc.Qj_0 / (twl * s0 * Ck * dx))))))

        else:
            Ck = 0
            X = 0.5

    if Ck > 0:
        Km = fmaxf(dt, dx/Ck)
    else:
        Km = dt

    D = (Km*(1.0 - X) + 0.5 * dt)              # --seconds

    qc.C1 = (Km*X + dt/2)/D
    qc.C2 = (dt/2 - Km*X)/D
    qc.C3 = (Km*(1-X)-dt/2)/D
    qc.C4 = (ql*dt)/D


    cdef float C_dot_Q

    if WP + WPC > 0:
        C_dot_Q = (qc.C1 * qup) + (qc.C2 * quc) + (qc.C3 * qdp) + qc.C4
        qc.Qj_0 = C_dot_Q - ((1 / (((WP * n) + (WPC * ncc)) / (WP + WPC))) *
                                                              (AREA + AREAC) * (powf(R,(2./ 3.))) * sqr_s0)



cdef inline void secant_h(
    float z,
		float bw,
		float bfd,
		float twcc,
		float s0,
		float n,
		float ncc,
		float dt,
		float dx,
		float qup,
		float quc,
		float qdp,
		float ql,
		float h,
    QC* qc,
) nogil:
    cdef float twl, AREA, AREAC, WP, WPC, R, Ck, X, Km, D
    cdef float areasum

    twl = bw + 2.0 * z * h

    cdef float sqr_s0 = sqrtf(s0)
    cdef float sqr_1z2 = sqrtf(1 + powf(z,2))
    cdef float C_dot_Q = (qc.C1 * qup) + (qc.C2 * quc) + (qc.C3 * qdp) + qc.C4

    if h > bfd:
        AREA = (bw + bfd * z) * bfd
        AREAC = (twcc * (h-bfd))
        WP = (bw + 2.0 * bfd * sqr_1z2)
        WPC = twcc + (2.0*(h-bfd))
        R = (AREA + AREAC)/(WP + WPC)

        areasum = 1/(AREA+AREAC)
        Ck = fmaxf(0.0, ((sqr_s0 / n) * ((5. / 3.) * powf(R,(2. / 3.)) -
                                         ((2. / 3.) * powf(R,(5. / 3.)) * (
                                                     2.0 * sqr_1z2 / (bw + 2.0 * bfd * z)))) * AREA
                       + ((sqr_s0 / (ncc)) * (5. / 3.) * powf((h - bfd),(2. / 3.))) * AREAC) *areasum)

        X = fminf(0.5, fmaxf(0.25, 0.5 * (1 - (C_dot_Q / (2 * twcc * s0 * Ck * dx)))))

    else:
        AREA = (bw + h * z ) * h
        AREAC = 0
        WP = (bw + 2.0 * h * sqr_1z2)
        WPC = 0
        if WP > 0:
            R = AREA/WP
        else:
            R = 0

        if h > 0:
            Ck = fmaxf(0.0, (sqr_s0 / n) * ((5. / 3.) * powf(R,(2. / 3.)) -
                                            ((2. / 3.) * powf(R,(5./ 3.)) * (
                                                        2.0 * sqr_1z2 / (bw + 2.0 * h * z)))))
        else:
            Ck = 0

        if Ck > 0:
            X = fminf(0.5, fmaxf(0.25, 0.5 * (1 - (C_dot_Q / (2 * twl * s0 * Ck * dx)))))
        else:
            X = 0.5

    if Ck > 0:
        Km = fmaxf(dt, dx/Ck)
    else:
        Km = dt
    D = 1/(Km*(1-X) + 0.5*dt)

    qc.C1 = (Km*X + 0.5 *dt) * D
    qc.C2 = (0.5 * dt - Km*X) * D
    qc.C3 = (Km*(1.0-X)-0.5 * dt) * D
    qc.C4 = (ql*dt) * D


    cdef float t

    t = (qc.C1 * qup) + (qc.C2 * quc) + (qc.C3 * qdp)
    if qc.C4 < 0 and fabsf(qc.C4) > t:
        qc.C4 = -t

    if WP + WPC > 0:
        qc.Qj = (t + qc.C4) - (((WP + WPC)/(((WP*n)+(WPC*ncc)))) *
                    (AREA+AREAC) * (powf(R,(2./3.))) * sqr_s0)
