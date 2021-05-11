//#include "reach.h"

static const double TWOTHIRDS= 2.0/3.0;
static const double FIVETHIRDS = 5.0/3.0;


typedef struct QVD {
    double qdc;
    double velc;
    double depthc;
    double cn;
    double ck;
    double X;
} QVD_double;


typedef struct QHC {
    double h;
    double Q_mc;
    double Q_normal;
    double Q_j;
    double Xmin;
    double X;
    double ck;
    double cn;
    double C1;
    double C2;
    double C3;
    double C4;
} QHC;

typedef struct channel_properties {
    double bfd;
    double bw;
    double tw;
    double twcc;
    double z;
    double s0;
    double sqrt_s0;
    double sqrt_1z2;
    double n;
    double ncc;
} channel_properties;


typedef struct hydraulic_geometry {
    double twl;
    double R;
    double AREA;
    double AREAC;
    double WP;
    double WPC;
    double h_lt_bf;
    double h_gt_bf;
} hydraulic_geometry;




void compute_hydraulic_geometry(const double, const channel_properties*, hydraulic_geometry*);
void compute_mc_flow(const channel_properties*, 
                     const double, 
                     const double, 
                     const double,
                     const double,
                     const double,
                     const double,
                     hydraulic_geometry*,
                     QHC*);
void compute_celerity(
    const channel_properties*,
    const hydraulic_geometry*,
    QHC*);

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
);
