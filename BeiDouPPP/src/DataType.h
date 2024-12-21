#pragma once
#include <math.h>
#include <time.h>
#include <stdint.h>


#define PI          3.1415926535897932  /* pi */

#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define GM          3.986004418E14      /*  standard gravitational parameter of earth   */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */

#define  OMGe         7.292115E-5    /* rad of earth's rotation*/
#define AU          149597870691.0      /* 1 AU (m) */
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */


#define FREQ1_CMP   1.561098E9          /* BeiDou B1 frequency (Hz) */
#define FREQ3_CMP   1.26852E9           /* BeiDou B3 frequency (Hz) */
#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */

#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   61                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */

#define MAXSAT   150           /* max satellite number (1 to MAXSAT) */ 
#define NFREQ        2                   /* number of carrier frequencies */

#define VAR_POS     SQR(60)      /* init variance receiver position (m^2) */
#define VAR_CLK      SQR(60)      /* init variance receiver clock (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */

#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */
#define THRES_MW_JUMP 10

#define NMAX        10              /* order of polynomial interpolation */

#define MAXLEAPS    64                  /* max number of leap seconds table */

#define SQR(x)      ((x)*(x))

#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
/* coordinate rotation matrix ------------------------------------------------*/
#define Rx(t,X) do { \
    (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
    (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
    (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
    (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

////////////PANLIN
const double BDS_IGSO[3][10] = { {-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},
                               {-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},
                               {-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32} };

const double BDS_MEO[3][10] = { {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},
                              {-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},
                              {-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47} };

const double leaps[MAXLEAPS + 1][7] = { /* leap seconds (y,m,d,h,m,s,utc-gpst) */
    {2017,1,1,0,0,0,-18},
    {2015,7,1,0,0,0,-17},
    {2012,7,1,0,0,0,-16},
    {2009,1,1,0,0,0,-15},
    {2006,1,1,0,0,0,-14},
    {1999,1,1,0,0,0,-13},
    {1997,7,1,0,0,0,-12},
    {1996,1,1,0,0,0,-11},
    {1994,7,1,0,0,0,-10},
    {1993,7,1,0,0,0, -9},
    {1992,7,1,0,0,0, -8},
    {1991,1,1,0,0,0, -7},
    {1990,1,1,0,0,0, -6},
    {1988,1,1,0,0,0, -5},
    {1985,7,1,0,0,0, -4},
    {1983,7,1,0,0,0, -3},
    {1982,7,1,0,0,0, -2},
    {1981,7,1,0,0,0, -1},
    {0}
};

const double range[] = { 0.00,360.00,-90.00,90.00 };



/* time struct */
struct gtime_t
{
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
};

/* dcb data (satellite) struct */
struct dcb_t
{
    CString sat;
    double bias;    
};

struct gpst_t
{
    int week;               /* gps week*/
    double sec;            /* time of week*/
};

/* observation data record */
struct obsd_t
{
    bool rej;                   /* rejection flag(0:accpet,1:reject) */
    int prn;                    /* prn number of satellite */
    bool svh;                   /* sv health flag*/
    int EpochIndex;        /* index of epoch */
    CString sat;              /* satellite/receiver number */
    gtime_t time;       /* receiver sampling time (GPST) */
    double P[NFREQ]; /* observation data pseudorange (m) */
    double L[NFREQ]; /* observation data carrier-phase (cycle) */
};

/* observation data */
struct obs_t
{
    int*validIndex;     /* the index of valid data in all data */
    int n,nmax;         /* number of obervation data/allocated*/
    obsd_t* data;       /* observation data records */
};

/* BeiDou broadcast ephemeris type */
struct eph_t
{
    CString sat;            /* satellite number */
    //int flag;           /* GPS/QZS: L2 P data flag */
    /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    gtime_t toc, ttr; /* Toe,Toc,T_trans */
    double toe;         /* Toe (s) in week */
    int sva;               /* SV accuracy (URA index) */
    /* SV orbit parameters */
    double sqrtA, e, i0, OMG0, omg, M0, deln, OMGd, idot;
    double crc, crs, cuc, cus, cic, cis;        
    double fit;         /* fit interval (h) */
    double f0, f1, f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[6];      /* group delay parameters */
    /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
    /*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
};

/* precise eph type */
struct peph_t
{
    CString sat;            /* satellite number */
    double pos[3]; /* satellite position (ecef) (m|s) */
    float  std[3]; /* satellite position std (m|s) */
    double vel[3]; /* satellite velocity-rate (m/s) */
    float  vst[3]; /* satellite velocity-rate std (m/s) */
    float  cov[3]; /* satellite position covariance (m^2) */
    float  vco[3]; /* satellite velocity covariance (m^2) */
};

/* precise ephs type */
struct pephs_t 
{   
    int EpochIndex;        /* index of epoch */
    gtime_t time;       /* time (GPST) */
    peph_t* satstate;   /* satellite state parameters */
    int SatCount;   /* satellite counts */
} ;

/* precise clock type */
struct pclk_t {
    CString sat;            /* satellite number */
    double clk; /* satellite clock (s) */
    float std; /* satellite clock std (s) */
};

/* precise clocks type */
struct pclks_t {
    gtime_t time;       /* time (GPST) */
    int n, nmax;         /* number of data/allocated */
    pclk_t* pclk;             /* precise clocks data */
};

/* antenna parameter type */
struct pcv_t
{
    int flag;                   /*1:satellite 0:receiver */
    CString sat;            /* satellite number */
    gtime_t ts, te;      /* valid time start and end */
    double off[NFREQ][3]; /* phase center offset (PCO) e/n/u or x/y/z (m) */
    double var[NFREQ][19]; /* phase center variation (PCV) (m) */
    /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
};

/* antenna parameters type */
struct pcvs_t
{
    int n, nmax;         /* number of data/allocated */
    pcv_t* pcv;         /* antenna parameters data */
};

/* earth rotation parameter data type */
struct erpd_t 
{       
    double mjd;         /* mjd (days) */
    double xp, yp;       /* pole offset (rad) */
    double xpr, ypr;     /* pole offset rate (rad/day) */
    double ut1_utc;     /* ut1-utc (s) */
    double lod;         /* length of day (s/day) */
} ;

/* navigation data type */
struct  nav_t
{
    int n, np, nc;         /* number of broadcast/precise ephemeris/clock */
    pclks_t* pclks;      /* precise clock */
    eph_t* eph;         /* BDS ephemeris */
    pephs_t* peph;   /* precise ephemeris */
    pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
    double utc_cmp[8];  /* BeiDou UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    dcb_t cbias[MAXSAT];/* satellite dcb (C2I-C6I) (m) */
    double rbias;/* receiver dcb (C2I-C6I) (m) */
};

/* solution type */
struct sol_t
{
    int chisqrTest;   /* Chi-square Test */
    gtime_t time;       /* time (GPST) */
    double rr[3];       /* position (m) */
    /* {x,y,z} or {e,n,u} */
    double  qr[6];       /* position variance/covariance (m^2) */
    /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
    /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    double dtr;      /* receiver clock bias to time systems (s) */
    double gdop, pdop;  /* gdop/pdop */
    int ns;         /* number of valid satellites */
};

/* solution buffer type */
struct solbuf_t
{
    int n, nmax;         /* number of solution/max number of buffer */
    int start, end;      /* start/end index */
    gtime_t time;       /* current solution time */
    sol_t* data;        /* solution data */
    double rb[3];       /* reference position {x,y,z} (ecef) (m) */
};

/* satellite status type */
struct ssat_t 
{        
    //uint8_t sys;        /* navigation system */
    //uint8_t vs;         /* valid satellite flag single */
    double azel[2];     /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ]; /* residuals of pseudorange (m) */
    double resc[NFREQ]; /* residuals of carrier-phase (m) */
    //uint8_t vsat[NFREQ]; /* valid satellite flag */
    //uint16_t snr[NFREQ]; /* signal strength (*SNR_UNIT dBHz) */
    //uint8_t fix[NFREQ]; /* ambiguity fix flag (1:fix,2:float,3:hold) */
    bool slip[NFREQ]; /* cycle-slip flag */
    //uint8_t half[NFREQ]; /* half-cycle valid flag */
    //uint8_t half[NFREQ]; /* half-cycle valid flag */
    //int lock[NFREQ];   /* lock counter of phase */
    //uint32_t outc[NFREQ]; /* obs outage counter of phase */
    //uint32_t slipc[NFREQ]; /* cycle-slip counter */
    //uint32_t rejc[NFREQ]; /* reject counter */
    double gf[NFREQ - 1]; /* geometry-free phase (m) */
    double mw[NFREQ - 1]; /* MW-LC (m) */
   double phw;         /* phase windup (cycle) */
    //gtime_t pt[2][NFREQ]; /* previous carrier-phase time */
    //double ph[2][NFREQ]; /* previous carrier-phase observable (cycle) */
} ;

/* rtk control/result type */
struct rtk_t 
{   
    sol_t  sol;         /* PPP solution */
    double tt;          /* time difference between current and previous (s) */
    int nx, na;          /* number of float states/fixed states */
    CMatrix x, P;    /* float states and their covariance */
    CMatrix xa, Pa;     /* fixed states and their covariance */
    ssat_t ssat[NSATCMP];   /* slip state of satellite */
} ;
