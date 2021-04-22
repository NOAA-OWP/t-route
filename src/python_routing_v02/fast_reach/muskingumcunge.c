#include <math.h>
#include <float.h>
#include "muskingumcunge.h"


void compute_hydraulic_geometry(
    const double h,
    const channel_properties *chan,
    hydraulic_geometry *hg) {

    hg->twl = chan->bw + 2*chan->z*h;

    hg->h_gt_bf = fmax(h - chan->bfd, 0);
    hg->h_lt_bf = fmin(chan->bfd, h);

    hg->AREA = (chan->bw + hg->h_lt_bf * chan->z) * hg->h_lt_bf;

    hg->WP = (chan->bw + 2 * hg->h_lt_bf * chan->sqrt_1z2);

    hg->AREAC = (chan->twcc * hg->h_gt_bf);

    if(hg->h_gt_bf > 0) {
        hg->WPC = chan->twcc + (2 * hg->h_gt_bf);
    } else {
        hg->WPC = 0;
    }

    hg->R = (hg->AREA + hg->AREAC)/(hg->WP + hg->WPC);
}


void compute_mc_flow(
    const channel_properties *chan,
    const double dt,
    const double dx,
    const double qup,
    const double quc,
    const double qdp,
    const double ql,
    hydraulic_geometry *hg,
    QHC *qc) {

    compute_hydraulic_geometry(qc->h, chan, hg);
    compute_celerity(chan, hg, qc);
    qc->cn = qc->ck * (dt/dx);

    double tw;
    if(qc->h > chan->bfd) {
        tw = chan->twcc;
    } else {
        tw = hg->twl;
    }

    const double Xmin = qc->Xmin;
    const double Xtmp = 0.5 * (1 - (qc->Q_j / (2 * tw * chan->s0 * qc->ck * dx)));
    if(Xtmp <= Xmin) {
       qc->X = Xmin;
    } else if((Xtmp > Xmin) && (Xtmp < 0.5)) {
       qc->X = Xtmp;
    } else {
       qc->X = 0.5;
    }

    const double Km = dx/qc->ck;
    if((qc->ck > 0) && (dt < Km)) {
        const double dt2 = dt/2.0;
        const double tmp1 = (dt - 2.0 * Km * qc->X);
        const double tmp2 = Km * (1 - qc->X);
        const double tmp3 = tmp1 + 2.0 * Km;
        qc->C1 = (dt + 2 * Km * qc->X) / tmp3;
        qc->C2 = tmp1 / tmp3;
        qc->C3 = (tmp2 - dt2) / (tmp2 + dt2);
        qc->C4 = (2 * ql * dt) / tmp3;
    } else {
        const double D = 1.5 - qc->X;
        qc->C1 = (qc->X + 0.5)/D;
        qc->C2 = (0.5 - qc->X)/D;
        qc->C3 = qc->C2;
        qc->C4 = ql/D;
    }

    const double t = (qc->C1 * qup) + (qc->C2 * quc) + (qc->C3 * qdp);
    /* qc.C4 cannot be more negative than the sum of other terms */
    qc->C4 = fmax(-t, qc->C4);

    if(hg->WP + hg->WPC > 0) {
        qc->Q_mc = t + qc->C4;
        qc->Q_normal = (
            (1.0 / (((hg->WP * chan->n) + (hg->WPC * chan->ncc)) / (hg->WP + hg->WPC)))
            * (hg->AREA + hg->AREAC)
            * pow(hg->R, TWOTHIRDS)
            * chan->sqrt_s0
        );
        qc->Q_j = qc->Q_mc - qc->Q_normal;
    } else {
        qc->Q_j = 0.0;
    }
}


void compute_celerity(
    const channel_properties *chan,
    const hydraulic_geometry *hg,
    QHC *qc) {

    if (qc->h > chan->bfd) {
        qc->ck = fmax(
            0.0,
            (
                (chan->sqrt_s0 / chan->n)
                * (
                    FIVETHIRDS * pow(hg->R, TWOTHIRDS)
                    - (
                        TWOTHIRDS * pow(hg->R, FIVETHIRDS)
                        * (
                            2.0
                            * chan->sqrt_1z2
                            / (chan->bw + 2.0 * chan->bfd * chan->z)
                        )
                    )
                )
                * hg->AREA
                + (
                    (chan->sqrt_s0 / chan->ncc)
                    * FIVETHIRDS * pow(qc->h - chan->bfd, TWOTHIRDS)
                )
                * hg->AREAC
            )
            / (hg->AREA + hg->AREAC)
        );
    } else if(qc->h > 0.0) {
        qc->ck = fmax(
            0.0,
            (chan->sqrt_s0 / chan->n)
            * (
                FIVETHIRDS * pow(hg->R, TWOTHIRDS)
                - (
                    TWOTHIRDS * pow(hg->R, FIVETHIRDS)
                    * (
                        2.0
                        * chan->sqrt_1z2
                        / (chan->bw + 2.0 * qc->h * chan->z)
                    )
                )
            )
        );
    } else {
        qc->ck = 0.0;
    }
}

void muskingum_cunge(
    double dt,
    double qup,
    double quc,
    double qdp,
    double ql,
    double dx,
    double bw,
    double tw,
    double twcc,
    double n,
    double ncc,
    double cs,
    double s0,
    double velp,
    double depthp,
    QVD_double *rv
) {

    int maxiter = 100;
    double mindepth = 0.01;
    double aerror = 0.01;
    double rerror = 1.0;

    int it = 0;
    int tries = 0;

    double depthc = fmax(0, depthp);

    channel_properties chan_struct;
    channel_properties *chan = &chan_struct;
    chan_struct.n = n;
    chan_struct.ncc = ncc;
    chan_struct.bw = bw;
    chan_struct.tw = tw;
    chan_struct.twcc = twcc;
    chan_struct.s0 = s0;
    chan_struct.sqrt_s0 = sqrt(s0);

    if(cs == 0) {
        chan_struct.z = 1.0;
        chan_struct.z = sqrt(2);
    } else {
        chan_struct.z = 1.0/cs;
        chan_struct.sqrt_1z2 = sqrt(1.0 + chan_struct.z * chan_struct.z);
    }

    if(chan_struct.bw > chan_struct.tw) {
        chan_struct.bfd = chan_struct.bw * (1/0.00001);
    } else if(chan_struct.bw == chan_struct.tw) {
        chan_struct.bfd = 0.5 * (chan_struct.bw/chan_struct.z);
    } else {
        chan_struct.bfd = 0.5 * (chan_struct.tw - chan_struct.bw)/chan_struct.z;
    }

    QHC qc_struct_left;
    QHC *qc_left = &qc_struct_left;
    qc_struct_left.h = depthp * TWOTHIRDS;
    qc_struct_left.Q_mc = 0.0;
    qc_struct_left.Q_normal = 0.0;
    qc_struct_left.Q_j = 0.0;
    qc_struct_left.Xmin = 0.0;
    qc_struct_left.X = 0.0;
    qc_struct_left.ck = 0.0;
    qc_struct_left.cn = 0.0;
    qc_struct_left.C1 = 0.0;
    qc_struct_left.C2 = 0.0;
    qc_struct_left.C3 = 0.0;
    qc_struct_left.C4 = 0.0;

    QHC qc_struct_right;
    QHC *qc_right = &qc_struct_right;
    qc_struct_right.h = (depthp * (4.0/3.0)) + mindepth;
    qc_struct_right.Q_mc = 0.0;
    qc_struct_right.Q_normal = 0.0;
    qc_struct_right.Q_j = 0.0;
    qc_struct_right.Xmin = 0.25;
    qc_struct_right.X = 0.0;
    qc_struct_right.ck = 0.0;
    qc_struct_right.cn = 0.0;
    qc_struct_right.C1 = 0.0;
    qc_struct_right.C2 = 0.0;
    qc_struct_right.C3 = 0.0;
    qc_struct_right.C4 = 0.0;

    hydraulic_geometry hg_struct;
    hydraulic_geometry *hg = &hg_struct;

    if (ql > 0 || quc > 0 || qup > 0 || qdp > 0) {
        double h, h_1, h_0;

        while ((rerror > 0.01) && (aerror >= mindepth) && (it <= maxiter)) {
            /* compute secant_h0 and secant_h */
            compute_mc_flow(chan, dt, dx, qup, quc, qdp, ql, hg, qc_left);
            qc_right->Q_j = qc_left->Q_mc;
            compute_mc_flow(chan, dt, dx, qup, quc, qdp, ql, hg, qc_right);

            h_0 = qc_left->h;
            h = qc_right->h;

            if(qc_left->Q_j - qc_right->Q_j > DBL_EPSILON) {
                h_1 = h - (qc_right->Q_j * (h_0 - h)) / (qc_left->Q_j - qc_right->Q_j);
                if(h_1 < 0) {
                    h_1 = h;
                }
            } else {
                h_1 = h;
            }

            if(h > 0) {
                rerror = fabs((h_1 - h) * (1/h));
                aerror = fabs(h_1 - h);
            } else {
                rerror = 0;
                aerror = 0.9;
            }

            it++;
            qc_left->h = fmax(0, h);
            qc_right->h = fmax(0, h_1);

            if (h < mindepth) {
                if (it >= maxiter) {
                    tries++;
                    if (tries <= 4) {
                        h = h * (4.0/3.0);
                        h_0 = h_0 * TWOTHIRDS;
                        maxiter += 25;
                        qc_left->Q_j = 0;
                        continue;
                    } else {
                        break;
                    }
                }
            }
        }

        h = qc_right->h;
        const double C1 = qc_right->C1;
        const double C2 = qc_right->C2;
        const double C3 = qc_right->C3;
        const double C4 = qc_right->C4;
        const double C_pdot_Q = (C1 * qup) + (C2 * quc) + (C3 * qdp);
        const double qdc = C_pdot_Q + C4;
        if ((qdc < 0) && (C4 < 0) && (fabs(C4) > C_pdot_Q)) {
            rv->qdc = 0;
        } else {
            rv->qdc = fmax((C1 * qup) + (C2 * quc) + C4, (C1*qup) + (C3 * qdp) + C4);
        }

        const double twl = bw + (2 * chan_struct.z * h);
        const double R = (0.5 * h * (bw+twl)) / (bw + 2.0 * sqrt(pow(((twl - bw) * 0.5), 2) + (h * h)));
        rv->velc = (1.0/n) * (pow(R, TWOTHIRDS)) * chan_struct.sqrt_s0;
        rv->depthc = qc_right->h;

        compute_celerity(chan, hg, qc_right);
        rv->ck = qc_right->ck;
        rv->cn = qc_right->ck * (dt / dx);
        rv->X = 0.0;
    } else {
        rv->qdc = 0.0;
        rv->depthc = 0.0;
        rv->velc = 0.0;
        rv->cn = 0.0;
        rv->ck = 0.0;
        rv->X = 0.0;
    }

}
