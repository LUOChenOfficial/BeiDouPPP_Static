#pragma once
#include "..//src/BdSPP.h"
#include "geoid.h"

class CBdPPP:public CBdSPP
{
protected:
	double odisp[6 * 11]; /* ocean tide loading parameters */
	erpd_t erp;	/* earth rotation parameter data */
	double prn[6];      /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv [5] pos */
	double thresslip;	/* threshold of slip detection */
	double maxinno;     /* reject threshold of innovation (m) */
	double antdel[3];  /* antenna delta {rov_e,rov_n,rov_u} */
	pcv_t pcvr;		/* pcv data of receiver */
	double eratio[NFREQ]; /* code/phase error ratio */
	double err[5];      /* measurement error factor */
	/* [0]:reserved */
	/* [1-3]:error factor a/b/c of phase (m) */
	/* [4]:doppler frequency (hz) */
public:
	CBdPPP();
	~CBdPPP();
public:
	/* READ DATA FUNCTION */
	/* read precise eph data	 */
	void ReadPreEphData(CString Sp3File);
	/* read precise clk data */
	void ReadPreClkData(CString ClkFile);
	/* read pcv and pco data */
	void ReadPcvData(CString AtxFile);
	/* read erp data */
	void ReadErpData(CString ErpFile);
	/* read otl data */
	void ReadOtlData(CString BlqFile);
	/* read dcb data */
	void ReadDcbData(CString BsxFile);

	/* PROCESS FUNCTION */
	/* initialize rtk control */
	void rtkinit(rtk_t& rtk);
	/* update ppp states */
	void udstate_ppp(obs_t obs);
	/* calculate satellite position of every epoch by PRECISE EPHEMERIS*/
	CMatrix CalPreSatPosClkEachEpoch(int epoIndex);
	/* tidal displacement */
	void tidedisp(gtime_t tutc, const double* rr, const double* odisp, double* dr);
	/* form design matrix and constant matrix (H&R&v) of KalmanFilter for PPP*/
	int FormKFDesConsMatEachEpoch(int post,int epoIndex, obs_t obs, CMatrix SatPosClkEpoch, CMatrix& x0, CMatrix& H, CMatrix& R, CMatrix& v,  double* dr);

	/* MAIN PROCESS */
	/* precise point position(ppp) */
	void pppos(CString testMat);

	/* Output PROCESS */
	void OutputSol(CString SolFile);
protected:
	/* TOOL FUNCTION */
	/* update postion state of ppp */
	void udpos_ppp();
	/* update clock state of ppp */
	void udclk_ppp();
	/* update tropospheric parameters of ppp*/
	void udtrop_ppp();
	/* update phase-bias of ppp*/
	void udbias_ppp(obs_t obs);
	/* tropospheric delay correction */
	double sbstropcorr(gtime_t time, const double* blh, const double* azel, double* var);
	/* time to day of year */
	double time2doy(gtime_t t);
	/* convert calendar day/time to time */
	gtime_t epoch2time(const double* ep);
	/* utc to gpstime */
	gtime_t utc2gpst(gtime_t t);
	/* gpstime to utc */
	gtime_t gpst2utc(gtime_t t);
	/* utc to gmst */
	double utc2gmst(gtime_t t, double ut1_utc);
	/* time to day and sec */
	double time2sec(gtime_t time, gtime_t* day);
	/* get meterological parameters */
	void getmet(double lat, double* met);
	/* detect cycle slip by geometry free phase jump */
	void detslp_gf(obs_t obs);
	/* detect slip by Melbourne-Wubbena linear combination jump */
	void detslp_mw(obs_t obs);
	/* geometry-free phase measurement */
	double gfmeas(obsd_t obs);
	/* Melbourne-Wubbena linear combination */
	double mwmeas(obsd_t obs);
	/* antenna corrected measurements */
	void corr_meas(obsd_t obs, const double* azel, const double* dantr, const double* dants, double phw, double* L, double* P, double* Lc, double* Pc);
	/* polynomial interpolation by Neville's algorithm */
	double interppol(const double* x, double* y, int n);
	/* get sat pos by precise eph */
	int pephpos(double* rs, gtime_t time, int epoIndex, int satIndex);
	/* get sat clk by precise eph */
	int pephclk(double& dts, gtime_t time,int epoIndex,int satIndex);
	/* satellite antenna phase center offset */
	void satantoff(gtime_t time, const double* rs, CString sat, double* dant);
	/* sun and moon position */
	void sunmoonpos(gtime_t tutc, const double* erpv, double* rsun, double* rmoon, double* gmst);
	/* sun and moon position in eci */
	void sunmoonpos_eci(gtime_t tut, double* rsun, double* rmoon);
	/* eci to ecef transformation matrix */
	void eci2ecef(gtime_t tutc, const double* erpv, double* U, double* gmst);
	/* astronomical arguments: f={l,l',F,D,OMG} (rad) */
	void ast_args(double t, double* f);
	/* iau 1980 nutation */
	void nut_iau1980(double t, const double* f, double* dpsi, double* deps);
	/* displacement by solid earth tide (ref [2] 7) */
	void tide_solid(const double* rsun, const double* rmoon, const double* pos, const double* E, double gmst, double* dr);
	/* displacement by ocean tide loading (ref [2] 7) */
	void tide_oload(gtime_t tut, const double* odisp, double* denu);
	/* displacement by pole tide (ref [7] eq.7.26) */
	void tide_pole(gtime_t tut, const double* pos, const double* erpv, double* denu);
	/* solar/lunar tides (ref [2] 7) */
	void tide_pl(const double* eu, const double* rp, double GMp, const double* pos, double* dr);
	/* iers mean pole (ref [7] eq.7.25) */
	void iers_mean_pole(gtime_t tut, double* xp_bar, double* yp_bar);
	/* tropospheric model */
	void model_trop(gtime_t time, const double* pos, const double* azel, const CMatrix x, double* dtdx, double* dtrp, double* var);
	/* satellite antenna phase center variation */
	void satantpcv(const double* rs, const double* rr, const pcv_t* pcv, double* dant);
	/* phase windup model */
	int model_phw(gtime_t time, CString sat, const double* rs, const double* rr, double* phw);
	/* precise tropospheric model */
	double trop_model_prec(gtime_t time, const double* pos, const double* azel, const double* x, double* dtdx, double* var);
	/* satellite antenna model */
	void antmodel_s(const pcv_t* pcv, double nadir, double* dant);
	/* receiver antenna model */
	void antmodel_r(const pcv_t* pcv, const double* del, const double* azel, double* dant);
	/* satellite attitude model */
	int sat_yaw(gtime_t time, CString sat, const double* rs, double* exs, double* eys);
	/* troposphere mapping function */
	double tropmapf(gtime_t time, const double pos[], const double azel[], double* mapfw);
	/* interpolate antenna phase center variation */
	double interpvar(double ang, const double* var);
	/* troposphere mapping function GMF */
	double gmf(gtime_t time, const double pos[], const double azel[], double* gmfw);
	/* yaw-angle of satellite */
	int yaw_angle(double beta, double mu, double* yaw);
	/* geoid height */
	double geoidh(const double* pos);
	/* embedded geoid model */
	double geoidh_emb(const double* pos);
	/* Legendre Spherical harmonics function */
	void Pnm(int n, int m, double x, double* P);
	double mapf(double el, double a, double b, double c);
	/* bilinear interpolation */
	double interpb(const double* y, double a, double b);
	/* outer product of 3d vectors */
	void cross3(const double* a, const double* b, double* c);
	/* nominal yaw-angle */
	double yaw_nominal(double beta, double mu);
	/* measurement error variance */
	double CBdPPP::varerr(CString sat, double el, int freq, int type, const double* pos);

	/* sqrt of var */
	double sqrtvar(double var);
	/* multiply matrix */
	void matmul(const char* tr, int n, int k, int m, double alpha, const double* A, const double* B, double beta, double* C);
};

