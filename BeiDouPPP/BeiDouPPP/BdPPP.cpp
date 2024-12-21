#include "pch.h"
#include "BdPPP.h"
CBdPPP::CBdPPP()
{
	erp = { 0 };
	prn[0] = 1e-4;
	prn[1] = 1e-2;
	prn[2] = 1e-4;
	prn[3] = 10;
	prn[4] = 10;
	prn[5] = 0;
	thresslip = 0.05;
	antdel[0] = 0;
	antdel[1] = 0;
	antdel[2] = 0;
	eratio[0] = eratio[1] = 100;
	err[0] = 0;
	err[1] = 0.003;
	err[2] = 0.003;
	err[3] = 0;
	err[4] = 10;
	navs.rbias = 42.594 * 1E-9 * CLIGHT;////////////read from dcb(station:JFNG)
	maxinno = 30;/////m
}

CBdPPP::~CBdPPP()
{
}

void CBdPPP::ReadPreEphData(CString Sp3File)
{
	CStdioFile sf;
	if (!sf.Open(Sp3File, CFile::modeRead)) return;
	CString strTmp;
	sf.ReadString(strTmp);
	int PreEphEpochCount = _ttoi(strTmp.Mid(37, 2));
	navs.peph = new pephs_t[PreEphEpochCount];
	navs.np = PreEphEpochCount;
	sf.ReadString(strTmp);
	sf.ReadString(strTmp);
	int SatCountTotal= _ttoi(strTmp.Mid(3, 3));

	/* read until data body*/
	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Right(2) != "/*");

	/* read body */
	for (int epoIndex = 0; epoIndex < PreEphEpochCount; epoIndex++)
	{
		sf.ReadString(strTmp);
		CString* strTime = new CString[6];
		gtime_t t;
		strTime[0] = strTmp.Mid(3, 4);
		strTime[1] = strTmp.Mid(9, 1);
		strTime[2] = strTmp.Mid(12, 1);
		strTime[3] = strTmp.Mid(14, 2).Trim();
		strTime[4] = strTmp.Mid(17, 2).Trim();
		strTime[5] = strTmp.Mid(21, 10);
		str2time(strTime, t);
		navs.peph[epoIndex].time = t;
		delete[]strTime;
		strTime = NULL;
		peph_t* pephTmp = new peph_t[SatCountTotal];
		int BeiDouSatCount = 0;
		for (int i = 0; i < SatCountTotal; i++)
		{
			sf.ReadString(strTmp);
			if (strTmp.Left(2) == "PC")BeiDouSatCount++;
			pephTmp[i].sat = strTmp.Mid(1, 3);
			pephTmp[i].pos[0] = _tstof(strTmp.Mid(5, 13).Trim()) * 1000;
			pephTmp[i].pos[1] = _tstof(strTmp.Mid(19, 13).Trim()) * 1000;
			pephTmp[i].pos[2] = _tstof(strTmp.Mid(33, 13).Trim()) * 1000;
		}
		
		/* extract BeiDou precise eph data */
		navs.peph[epoIndex].satstate = new peph_t[BeiDouSatCount];
		navs.peph[epoIndex].SatCount = BeiDouSatCount;
		for (int i = 0, j = 0; j < BeiDouSatCount; i++)
		{
			if (pephTmp[i].sat.Left(1) == 'C')
			{
				navs.peph[epoIndex].satstate[j].sat = pephTmp[i].sat;
				navs.peph[epoIndex].satstate[j].pos[0] = pephTmp[i].pos[0];
				navs.peph[epoIndex].satstate[j].pos[1] = pephTmp[i].pos[1];
				navs.peph[epoIndex].satstate[j].pos[2] = pephTmp[i].pos[2];
				j++;
			}
		}
		delete[]pephTmp;
		pephTmp = NULL;
	}
	sf.Close();

}

void CBdPPP::ReadPreClkData(CString ClkFile)
{
	CStdioFile sf;
	if (!sf.Open(ClkFile, CFile::modeRead)) return;
	CString strTmp;
	sf.ReadString(strTmp);
	sf.ReadString(strTmp);

	/* read time interval */
	double ti = _tstof(strTmp.Mid(4, 4).Trim());
	
	/* read until sat count */
	while (strTmp.Right(14) != "# OF SOLN SATS")
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	}
	int SatCountTotal = _ttoi(strTmp.Mid(0, 3));
	int ClkEpochCount = 24 * 3600 / ti ;
	navs.pclks = new pclks_t[ClkEpochCount];
	navs.nc = ClkEpochCount;

	/* read until data body */
	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Right(13) != "END OF HEADER");

	/* read body */
	//int BeiDouClkCount = 0;
	for (int epoIndex = 0; epoIndex < ClkEpochCount; epoIndex++) 
	{
		sf.ReadString(strTmp);

		if (epoIndex == 120)
		{
			double a = 1;
		}

		CString* strTime = new CString[6];
		gtime_t t;
		strTime[0] = strTmp.Mid(8, 4);
		strTime[1] = strTmp.Mid(14, 1);
		strTime[2] = strTmp.Mid(17, 1);
		strTime[3] = strTmp.Mid(19, 2).Trim();
		strTime[4] = strTmp.Mid(22, 2).Trim();
		strTime[5] = strTmp.Mid(25, 9).Trim();
		str2time(strTime, t);
		navs.pclks[epoIndex].time = t;
		delete[]strTime;
		strTime = NULL;

		navs.pclks[epoIndex].nmax = SatCountTotal;
		navs.pclks[epoIndex].pclk = new pclk_t[SatCountTotal];
		navs.pclks[epoIndex].pclk[0].sat = strTmp.Mid(3, 3);
		navs.pclks[epoIndex].pclk[0].clk = _tstof(strTmp.Mid(40, 19).Trim());
		for (int satInd = 1; satInd < SatCountTotal; satInd++)
		{
			sf.ReadString(strTmp);
			navs.pclks[epoIndex].pclk[satInd].sat = strTmp.Mid(3, 3);
			navs.pclks[epoIndex].pclk[satInd].clk = _tstof(strTmp.Mid(40, 19).Trim());
		}
	}
}

void CBdPPP::ReadPcvData(CString AtxFile)
{
	CStdioFile sf;
	if (!sf.Open(AtxFile, CFile::modeRead)) return;
	CString strTmp;

	int BeiDouCount = 69;
	int i;
	/* read pco/pcv of satellite */
	for(i=0;i<BeiDouCount;i++)
	{
		do
		{
			sf.ReadString(strTmp);
			strTmp.Trim();
		} while (strTmp.Left(6) != "BEIDOU" || strTmp.Right(9) != "SERIAL NO");

		navs.pcvs[i].sat = strTmp.Mid(20, 3);
		navs.pcvs[i].flag = 1;

		sf.ReadString(strTmp);
		sf.ReadString(strTmp);
		sf.ReadString(strTmp);
		sf.ReadString(strTmp);

		sf.ReadString(strTmp);
		CString* strTimeStart = new CString[6];
		gtime_t t;
		strTimeStart[0] = strTmp.Mid(2, 4);
		strTimeStart[1] = strTmp.Mid(10, 2).Trim();
		strTimeStart[2] = strTmp.Mid(16, 2).Trim();
		strTimeStart[3] = strTmp.Mid(22, 2).Trim();
		strTimeStart[4] = strTmp.Mid(28, 2).Trim();
		strTimeStart[5] = strTmp.Mid(33, 9).Trim();
		str2time(strTimeStart, t);
		navs.pcvs[i].ts = t;
	
		/* default end time (i.e. no " valid until ") */
		CString strTimeMax[] = { _T("2030"),_T("1"),_T("1"),_T("1"),_T("0"),_T("0") };
		str2time(strTimeMax, t);
		navs.pcvs[i].te = t;

		sf.ReadString(strTmp);
		if (strTmp.Trim().Right(11) == "VALID UNTIL")
		{
		CString* strTimeEnd = new CString[6];
		strTimeEnd[0] = strTmp.Mid(0, 4);
		strTimeEnd[1] = strTmp.Mid(8, 2).Trim();
		strTimeEnd[2] = strTmp.Mid(14, 2).Trim();
		strTimeEnd[3] = strTmp.Mid(20, 2).Trim();
		strTimeEnd[4] = strTmp.Mid(26, 2).Trim();
		strTimeEnd[5] = strTmp.Mid(31, 9).Trim();
		str2time(strTimeEnd, t);
		navs.pcvs[i].te = t;
		}


		do
		{
			sf.ReadString(strTmp);
			strTmp.Trim();
		} while (strTmp.Left(3) != "C02" || strTmp.Right(18) != "START OF FREQUENCY");

		sf.ReadString(strTmp);
		navs.pcvs[i].off[0][1] = _tstof(strTmp.Mid(14, 6).Trim()) * 1e-3;
		navs.pcvs[i].off[0][0] = _tstof(strTmp.Mid(3, 7).Trim()) * 1e-3;
		navs.pcvs[i].off[0][2] = _tstof(strTmp.Mid(23, 7).Trim()) * 1e-3;
		for (int ind = 0; ind < 19; ind++)navs.pcvs[0].var[0][ind] = 0 * 1e-3;

		do
		{
			sf.ReadString(strTmp);
			strTmp.Trim();
		} while (strTmp.Left(3) != "C06" || strTmp.Right(18) != "START OF FREQUENCY");

		sf.ReadString(strTmp);
		navs.pcvs[i].off[1][1] = _tstof(strTmp.Mid(14, 6).Trim()) * 1e-3;
		navs.pcvs[i].off[1][0] = _tstof(strTmp.Mid(3, 7).Trim()) * 1e-3;
		navs.pcvs[i].off[1][2] = _tstof(strTmp.Mid(23, 7).Trim()) * 1e-3;
		for (int ind = 0; ind < 19; ind++)navs.pcvs[0].var[0][ind] = 0 * 1e-3;

		do
		{
			sf.ReadString(strTmp);
			strTmp.Trim();
		} while (strTmp.Right(14) != "END OF ANTENNA");
	 }
	/* read pcv/pco of receiver */
	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Left(20) != "TRM59800.00     NONE");

	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Left(3) != "C02" || strTmp.Right(18) != "START OF FREQUENCY");
	sf.ReadString(strTmp);
	pcvr.off[0][0] = _tstof(strTmp.Mid(5, 5).Trim()) * 1e-3;
	pcvr.off[0][1] = _tstof(strTmp.Mid(15, 5).Trim()) * 1e-3;
	pcvr.off[0][2] = _tstof(strTmp.Mid(23, 6).Trim()) * 1e-3;
	sf.ReadString(strTmp);
	pcvr.var[0][0] = 0;
	for (int ind = 1; ind < 19; ind++)pcvr.var[0][ind] = _tstof(strTmp.Mid(19 + (ind - 1) * 8, 5).Trim()) * 1e-3;

	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Left(3) != "C06" || strTmp.Right(18) != "START OF FREQUENCY");
	sf.ReadString(strTmp);
	pcvr.off[1][0] = _tstof(strTmp.Mid(5, 5).Trim()) * 1e-3;
	pcvr.off[1][1] = _tstof(strTmp.Mid(15, 5).Trim()) * 1e-3;
	pcvr.off[1][2] = _tstof(strTmp.Mid(23, 6).Trim()) * 1e-3;
	sf.ReadString(strTmp);
	pcvr.var[1][0] = 0;
	for (int ind = 1; ind < 19; ind++)pcvr.var[1][ind] = _tstof(strTmp.Mid(19 + (ind - 1) * 8, 5).Trim()) * 1e-3;

	sf.Close();
}

void CBdPPP::ReadErpData(CString ErpFile)
{
	CStdioFile sf;
	if (!sf.Open(ErpFile, CFile::modeRead)) return;
	CString strTmp;

	sf.ReadString(strTmp);
	sf.ReadString(strTmp);
	sf.ReadString(strTmp);
	sf.ReadString(strTmp);
	sf.ReadString(strTmp);
	erp.mjd = _tstof(strTmp.Mid(0, 8));
	erp.xp = _tstof(strTmp.Mid(11, 5)) * 1e-6 / 3600.0 * PI / 180;
	erp.yp = _tstof(strTmp.Mid(18, 6)) * 1e-6 / 3600.0 * PI / 180;
	erp.xpr = _tstof(strTmp.Mid(81, 4)) * 1e-6 / 3600.0 * PI / 180;
	erp.ypr = _tstof(strTmp.Mid(88, 4)) * 1e-6 / 3600.0 * PI / 180;
	erp.ut1_utc = _tstof(strTmp.Mid(25, 8)) * 1e-7;
	erp.lod = _tstof(strTmp.Mid(35, 5)) * 1e-7;
}

void CBdPPP::ReadOtlData(CString BlqFile)
{
	CStdioFile sf;
	if (!sf.Open(BlqFile, CFile::modeRead)) return;
	CString strTmp;

	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Left(7) != "$$ JFNG");
	for (int i = 0; i < 6; i++)
	{
		sf.ReadString(strTmp);
		for (int j = 0; j < 11; j++)odisp[i + j * 6] = _tstof(strTmp.Mid(2 + 7 * j, 6).Trim());
	}
	int a = 1;
}

void CBdPPP::ReadDcbData(CString BsxFile)
{
	CStdioFile sf;
	if (!sf.Open(BsxFile, CFile::modeRead)) return;
	CString strTmp;

	do
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	} while (strTmp.Left(3) != "DSB" );

	int i = 0, count = 0;
	while (strTmp.Left(3) == "DSB")
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
		if (strTmp.Mid(10, 1) == "C" && strTmp.Mid(24, 8) == "C2I  C6I"&& strTmp.Mid(14, 4)=="    ")
		{
			navs.cbias[i].sat = strTmp.Mid(10, 3);
			navs.cbias[i].bias = _tstof(strTmp.Mid(82, 8).Trim()) * 1E-9 * CLIGHT;
			i++;
		}
		count++;
	 }
	count;
}



void CBdPPP::rtkinit(rtk_t& rtk)
{
	sol_t sol0 = { {0} };
	ssat_t ssat0 = { 0 };
	int i;

	rtk.sol = sol0;
	rtk.nx = 3 + 1 + 3 + NSATCMP;
	rtk.na = 3 + 1 + 3 + NSATCMP;
	rtk.tt = 0.0;
	rtk.x.SetSize(rtk.nx, 1);
	rtk.P.SetSize(rtk.nx, rtk.nx);
	rtk.xa.SetSize(rtk.na, 1);
	rtk.Pa.SetSize(rtk.na, rtk.na);
	for (i = 0; i < NSATCMP; i++) 
	{
		rtk.ssat[i] = ssat0;
	}
}

void CBdPPP::udstate_ppp( obs_t obs)
{
	/* update position */
	udpos_ppp();
	/* update clock */
	udclk_ppp();
	/* update tropospheric parameters */
	udtrop_ppp();
	/* update phase-bias */
	udbias_ppp(obs);
}

CMatrix CBdPPP::CalPreSatPosClkEachEpoch(int epoIndex)
{
	//////////最后一位（第8位）保存prn号[x/y/z+dt+vx/vy/vz+prn]
	CMatrix SatPosClkMax(obss[epoIndex].nmax, 8);
	CMatrix SatPosClkTmp(obss[epoIndex].n, 8);
	CMatrix SatPosClkAllBd(NSATCMP, 8);

	/* calculate detlta clock and position of every satellite */
	for (int satIndex = 0; satIndex < obss[epoIndex].nmax; satIndex++)
	{
		/* detlta clock */
		gtime_t ts, tr;
		double P, f0, f1, f2;
		double dt0, dt, dts;
		eph_t eph = seleph(obss[epoIndex].data[satIndex]);
		if (eph.e == 0)
		{
			obss[epoIndex].data[satIndex].rej = 1;
			SatPosClkMax(satIndex, 0) = 0;
			SatPosClkMax(satIndex, 1) = 0;
			SatPosClkMax(satIndex, 2) = 0;
			SatPosClkMax(satIndex, 3) = 0;
			SatPosClkMax(satIndex, 4) = 0;
			SatPosClkMax(satIndex, 5) = 0;
			SatPosClkMax(satIndex, 6) = 0;
			SatPosClkMax(satIndex, 7) = obss[epoIndex].data[satIndex].prn;
			continue;
		}
		double tk;
		gpst_t ts_gpst;

		tr = obss[epoIndex].data[satIndex].time;
		P = obss[epoIndex].data[satIndex].P[0];
		f0 = eph.f0;
		f1 = eph.f1;
		f2 = eph.f2;
		ts = timeadd(tr, -P / CLIGHT);
		dt0 = timediff(ts, eph.toc);
		double dtTmp = dt0;
		for (int i = 0; i < 2; i++)dt0 = dtTmp - (f0 + f1 * dt0 + f2 * dt0 * dt0);
		dt = f0 + f1 * dt0 + f2 * dt0 * dt0;

		/* sat postion */
		/* sat clock bias correction */
		ts = timeadd(ts, -dt);

		if (epoIndex == 60&& obss[epoIndex].data[satIndex].prn == 21)
		{
			double a = 1;
		}
		double  rs[6];
		double dts0 = 0, dts1 = 0, rs0[3], rs1[3];
		double tt = 1e-3;
		///////////////////////////////GET POS/////if peph of sat not found, continue to next
		if (!pephpos(rs0, ts, epoIndex, satIndex))
		{
			//obss[epoIndex].n = obss[epoIndex].n - 1;
			obss[epoIndex].data[satIndex].rej = 1;
			for (int i = 0; i < 7; i++)SatPosClkMax(satIndex, i) = 0;
			SatPosClkMax(satIndex, 7) = obss[epoIndex].data[satIndex].prn - 1;
			continue;
		};
		pephpos(rs1, timeadd(ts,tt), epoIndex, satIndex);
		///////////////////////////////GET CLK
		pephclk(dts0, ts, epoIndex, satIndex);
		pephclk(dts1, timeadd(ts, tt), epoIndex, satIndex);

		/* satellite antenna offset correction */
		CString sat = obss[epoIndex].data[satIndex].sat;
		double dant[3];

		satantoff(ts, rs0, sat, dant);

		/* difference sat pos to get velocity */
		for (int i = 0; i < 3; i++) 
		{
			rs[i] =rs0[i]+ dant[i];
			rs[i + 3] = (rs1[i] - rs0[i]) / tt;
		}

		/* relativistic effect correction */
		if (dts0 != 0.0) 
		{
			dts = dts0 - 2.0 * (rs[0] * rs[3] + rs[1] * rs[4] + rs[2] * rs[5]) / CLIGHT / CLIGHT;
		}

		SatPosClkMax(satIndex, 0) = rs[0];
		SatPosClkMax(satIndex, 1) = rs[1];
		SatPosClkMax(satIndex, 2) = rs[2];
		SatPosClkMax(satIndex, 3) = dts;
		SatPosClkMax(satIndex, 4) = rs[3];
		SatPosClkMax(satIndex, 5) = rs[4];
		SatPosClkMax(satIndex, 6) = rs[5];
		SatPosClkMax(satIndex, 7) = obss[epoIndex].data[satIndex].prn - 1;
	}		
	
	/* extract valid satellites */
	for (int i = 0; i < SatPosClkTmp.Row(); i++)
	{
		SatPosClkTmp(i, 0) = SatPosClkMax(obss[epoIndex].validIndex[i], 0);
		SatPosClkTmp(i, 1) = SatPosClkMax(obss[epoIndex].validIndex[i], 1);
		SatPosClkTmp(i, 2) = SatPosClkMax(obss[epoIndex].validIndex[i], 2);
		SatPosClkTmp(i, 3) = SatPosClkMax(obss[epoIndex].validIndex[i], 3);
		SatPosClkTmp(i, 4) = SatPosClkMax(obss[epoIndex].validIndex[i], 4);
		SatPosClkTmp(i, 5) = SatPosClkMax(obss[epoIndex].validIndex[i], 5);
		SatPosClkTmp(i, 6) = SatPosClkMax(obss[epoIndex].validIndex[i], 6);
		SatPosClkTmp(i, 7) = SatPosClkMax(obss[epoIndex].validIndex[i], 7);
	}

	CMatrix SatPosClkValid(obss[epoIndex].n, 8);
	for (int i = 0, j = 0; i < SatPosClkTmp.Row(); i++)
	{
		if (SatPosClkTmp(i, 0) == 0)continue;
			SatPosClkValid(j, 0) = SatPosClkTmp(i, 0);
			SatPosClkValid(j, 1) = SatPosClkTmp(i, 1);
			SatPosClkValid(j, 2) = SatPosClkTmp(i, 2);
			SatPosClkValid(j, 3) = SatPosClkTmp(i, 3);
			SatPosClkValid(j, 4) = SatPosClkTmp(i, 4);
			SatPosClkValid(j, 5) = SatPosClkTmp(i, 5);
			SatPosClkValid(j, 6) = SatPosClkTmp(i, 6);
			SatPosClkValid(j, 7) = SatPosClkTmp(i, 7);
			int ind = obss[epoIndex].validIndex[i];
			obss[epoIndex].validIndex[j] = ind;
			j++;
	}

	for (int i = 0; i < SatPosClkValid.Row(); i++)
	{
		int prnInd = SatPosClkValid(i, 7);
		SatPosClkAllBd(prnInd, 0) = SatPosClkValid(i, 0);
		SatPosClkAllBd(prnInd, 1) = SatPosClkValid(i, 1);
		SatPosClkAllBd(prnInd, 2) = SatPosClkValid(i, 2);
		SatPosClkAllBd(prnInd, 3) = SatPosClkValid(i, 3);
		SatPosClkAllBd(prnInd, 4) = SatPosClkValid(i, 4);
		SatPosClkAllBd(prnInd, 5) = SatPosClkValid(i, 5);
		SatPosClkAllBd(prnInd, 6) = SatPosClkValid(i, 6);
	}

	return SatPosClkAllBd;
}

void CBdPPP::tidedisp(gtime_t tutc, const double* rr, const double* odisp, double* dr)
{
	gtime_t tut;
	double B, L , E[9], drt[3], denu[3], rs[3], rm[3], gmst, erpv[5] = { 0 },pos[2];
	int i;

    //geterp(erp, utc2gpst(tutc), erpv);
	erpv[0] = erp.xp;
	erpv[1] = erp.yp;
	erpv[2] = erp.ut1_utc;
	erpv[3] = erp.lod;

	tut = timeadd(tutc, erpv[2]);

	dr[0] = dr[1] = dr[2] = 0.0;

	double norm = sqrt(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]);
	pos[0] = asin(rr[2] / norm);
	pos[1] = atan2(rr[1], rr[0]);
	B = pos[0];
	L = pos[1];

	E[0] = -sin(L);					  E[3] = cos(L);					E[6] = 0.0;
	E[1] = -sin(B) * cos(L);    E[4] = -sin(B) * sin(L);     E[7] = cos(B);
	E[2] = cos(B) * cos(L);    E[5] = cos(B) * sin(L);      E[8] = sin(B);

	 /* solid earth tides */
	{
		/* sun and moon position in ecef */
		sunmoonpos(tutc, erpv, rs, rm, &gmst);
		tide_solid(rs, rm, pos, E, gmst, drt);
		for (i = 0; i < 3; i++) dr[i] += drt[i];
	}
	 /* ocean tide loading */
	{
		tide_oload(tut, odisp, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i < 3; i++) dr[i] += drt[i];
	}
	  /* pole tide */
	{
		tide_pole(tut, pos, erpv, denu);
		matmul("TN", 3, 1, 3, 1.0, E, denu, 0.0, drt);
		for (i = 0; i < 3; i++) dr[i] += drt[i];
	}

}

int CBdPPP::FormKFDesConsMatEachEpoch(int post, int epoIndex, obs_t obs,CMatrix SatPosClkEpoch, CMatrix& x1, CMatrix& H, CMatrix& R, CMatrix& v, double* dr)
{
	double P[NFREQ], L[NFREQ];
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	double rr[3],pos[3];
	double obsi[NSATCMP * 2] = { 0 }, ve[NSATCMP * 2] = { 0 };
	int ne = 0, rej, maxobs;
	double vmax;
	int stat = 1;

	/* tidedisp correction */
	for (int i = 0; i < 3; i++)rr[i] = x1(i, 0) + dr[i];
	ecef2blh(rr, pos);

	/* form H&R&v for each satellite */
	for (int obsInd = 0; obsInd < obs.nmax; obsInd++)
	{
		int satInd = obs.data[obsInd].prn - 1;
		if (obs.data[obsInd].rej == 1 || SatPosClkEpoch(satInd, 0) == 0)
		{

			v(2 * obsInd, 0) = v(2 * obsInd + 1, 0) = 0;
			R(2 * obsInd, 2 * obsInd) = R(2 * obsInd + 1, 2 * obsInd + 1) = 0;
			for (int i = 0; i < H.Col(); i++)
			{
				H(2 * obsInd, i) = 0;
				H(2 * obsInd + 1, i) = 0;
			}
			continue;
		}
		/* search valid satellite in data index */
		CString sat = obs.data[obsInd].sat;
		double Lc, Pc;
		double dist, vec[3], SagEffectCor, azel[2];
		double rs[6], dts, cdtr = x1(3, 0);
		double dtrp = 0, varTrp = 0, dtdx[3];
		double bias = x1(satInd + 7, 0);
		pcv_t pcv_s = { 0 };
		for (int i = 0; i < 3; i++)rs[i] = SatPosClkEpoch(satInd, i);
		for (int i = 4; i < 7; i++)rs[i - 1] = SatPosClkEpoch(satInd, i);
		dts = SatPosClkEpoch(satInd, 3);

		/* tropospheric model */
		dist = geodist(rr, rs, vec, SagEffectCor);
		satazel(rr, vec, azel);
		if (fabs(azel[1]) < minElevAngle(RAD))
		{
			v(2 * obsInd, 0) = v(2 * obsInd + 1, 0) = 0;
			R(2 * obsInd, 2 * obsInd) = R(2 * obsInd + 1, 2 * obsInd + 1) = 0;
			for (int i = 0; i < H.Col(); i++)
			{
				H(2 * obsInd, i) = 0;
				H(2 * obsInd + 1, i) = 0;
			}
			continue;
		}
		dist += SagEffectCor;
		model_trop(obs.data[0].time, pos, azel, x1, dtdx, &dtrp, &varTrp);
		/* ionospheric model */
		////ILFC model (dion = 0)

		/* satellite and receiver antenna model */
		/* search for pcvs data */
		for (int i = 0; i < MAXSAT; i++)
		{
			double dts, dte;
			if (navs.pcvs[i].sat == sat && (timediff(obs.data[0].time, navs.pcvs[i].ts) > 0 && timediff(obs.data[0].time, navs.pcvs[i].te) < 0))
			{
				pcv_s = navs.pcvs[i];
				break;
			}
		}
		satantpcv(rs, rr, &pcv_s, dants);
		antmodel_r(&pcvr, antdel, azel, dantr);

		/* phase windup model */
		rtk.ssat[satInd].phw = 0;
		model_phw(rtk.sol.time, sat, rs, rr, &rtk.ssat[satInd].phw);

		/* corrected phase and code measurements */
		corr_meas(obs.data[obsInd], azel, dantr, dants, rtk.ssat[satInd].phw, L, P, &Lc, &Pc);

		/* form design matrix (H) */
		////////////Lc-H
		for (int i = 0; i < 3; i++)H(2 * obsInd, i) = -vec[i] / dist;
		H(2 * obsInd, 3) = 1.0;
		for (int i = 4; i < 7; i++)	H(2 * obsInd, i) = dtdx[i - 4];
		H(2 * obsInd, 7 + satInd) = 1.0;
		////////////Pc-H
		for (int i = 0; i < 3; i++)H(2 * obsInd + 1, i) = -vec[i] / dist;
		H(2 * obsInd + 1, 3) = 1.0;
		for (int i = 4; i < 7; i++)	H(2 * obsInd + 1, i) = dtdx[i - 4];
		/* form residual matrix (v) */
		v(2 * obsInd, 0) = Lc - (dist + cdtr - dts * CLIGHT + dtrp + bias);
		v(2 * obsInd + 1, 0) = Pc - (dist + cdtr - dts * CLIGHT + dtrp);
		/* form variance matrix (R) */
		R(2 * obsInd, 2 * obsInd) = varerr(sat, azel[1], 0, 0, pos) + varTrp;
		R(2 * obsInd + 1, 2 * obsInd + 1) = varerr(sat, azel[1], 0, 1, pos) + varTrp;

		/* reject satellite by pre-fit residuals */
		if (!post && (fabs(v(2 * obsInd, 0)) > maxinno * 2 || fabs(v(2 * obsInd + 1, 0)) > maxinno * 2))
		{
			v(2 * obsInd, 0) = v(2 * obsInd + 1, 0) = 0;
			R(2 * obsInd, 2 * obsInd) = R(2 * obsInd + 1, 2 * obsInd + 1) = 0;
			for (int i = 0; i < H.Col(); i++)
			{
				H(2 * obsInd, i) = 0;
				H(2 * obsInd + 1, i) = 0;
			}
			obs.data[obsInd].rej = 1;
			continue;
		}

		/* record large post-fit residuals *///////by Lc and Pc
		if (post) 
		{
			if (fabs(v(2 * obsInd, 0)) > sqrt(R(2 * obsInd, 2 * obsInd)) * THRES_REJECT)
			{
				obsi[ne] = 2 * obsInd; ve[ne] = v(2 * obsInd, 0); ne++;
			}
			if (fabs(v(2 * obsInd + 1, 0)) > sqrt(R(2 * obsInd + 1, 2 * obsInd + 1)) * THRES_REJECT)
			{
				obsi[ne] = 2 * obsInd + 1; ve[ne] = v(2 * obsInd + 1, 0); ne++;
			}
		}
	}
	/* reject satellite with large and max post-fit residual */
	if (post && ne > 0) 
	{
		vmax = ve[0]; maxobs = obsi[0]; rej = 0;
		for (int i = 1; i < ne; i++) 
		{
			if (fabs(vmax) >= fabs(ve[i])) continue;
			vmax = ve[i]; maxobs = obsi[i]; rej = i;
		}
		obs.data[maxobs / 2].rej = 1;
		ve[rej] = 0;
		stat = 0;
	}
	return post ? stat : 1;
}

void CBdPPP::OutputSol(CString SolFile)
{
	CStdioFile sf;
	if (!sf.Open(SolFile, CFile::modeCreate || CFile::modeNoTruncate)) return;
	//if (!sf.Open(SolFile, CFile::modeWrite)) return;
	CString strLine;
	// 查找从右向左第一个遇到的反斜杠 '\'
	//int pos1 = strObsFile.ReverseFind(_T('\\'));
	//CString strObsFileNotPath = strObsFile.Mid(pos1 + 1);
	//int pos2 = strNavFile.ReverseFind(_T('\\'));
	//CString strNavFileNotPath = strNavFile.Mid(pos2 + 1);

	/* write header*/

	strLine.Format(_T("%% obs file  : "));
	//strLine += strObsFileNotPath;
	strLine += _T("\r");
	sf.WriteString(strLine);
	strLine.Empty();

	strLine.Format(_T("%% nav file  : "));
	//strLine += strNavFileNotPath;
	strLine += _T("\r");
	sf.WriteString(strLine);
	strLine.Empty();

	strLine.Format(_T("%% elev mask : %.3f deg\r"), minElevAngle(DEG));
	sf.WriteString(strLine);
	strLine.Empty();
	sf.WriteString(_T("% ionos opt : broadcast\r"));
	sf.WriteString(_T("% tropo opt : saastamoinen\r"));
	sf.WriteString(_T("% ephemeris : broadcast\r"));
	sf.WriteString(_T("% navi sys  : beidou\r"));
	sf.WriteString(_T("%                       GPST      x-ecef(m)      y-ecef(m)      z-ecef(m)     clock bias  ns   sdx(m)   sdy(m)   sdz(m)  sdxy(m)  sdyz(m)  sdzx(m)   pdop   gdop chisqrTest\r"));

	/* write body*/
	for (int epoIndex = 0; epoIndex < EpochCount; epoIndex++)
	{
		CString strTmp;
		double ep[6];
		time2epoch(sol.data[epoIndex].time, ep);
		strLine.Format(_T("%4.0f\t%02.0f\t%02.0f\t%02.0f\t%02.0f\t%5.3f\t"), ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
		strTmp.Format(_T("%.4f\t%.4f\t%.4f\t%.4f\t  %d   "), sol.data[epoIndex].rr[0], sol.data[epoIndex].rr[1], sol.data[epoIndex].rr[2], sol.data[epoIndex].dtr*1E9, sol.data[epoIndex].ns);
		strLine += strTmp;
		strTmp.Format(_T("%.4f   %.4f\t%.4f\t%.4f   %.4f  %.4f    "), sqrt(sol.data[epoIndex].qr[0]), sqrt(sol.data[epoIndex].qr[1]), sqrt(sol.data[epoIndex].qr[2]), sqrtvar(sol.data[epoIndex].qr[3]), sqrtvar(sol.data[epoIndex].qr[4]), sqrtvar(sol.data[epoIndex].qr[5]));
		strLine += strTmp;
		strTmp.Format(_T("%.2f   %.2f\t%d\r"), sol.data[epoIndex].pdop, sol.data[epoIndex].gdop, sol.data[epoIndex].chisqrTest);
		strLine += strTmp;
		sf.WriteString(strLine);
		strTmp.Empty();
		strLine.Empty();
	}
	sf.Close();
}

void CBdPPP::pppos(CString testMat)
{
	rtkinit(rtk);
	gtime_t time;
	for (int epoIndex = 0; epoIndex < EpochCount; epoIndex++)
	{
		time = rtk.sol.time;
		pntpos(epoIndex);

		CMatrix SatPosClk = CalPreSatPosClkEachEpoch(epoIndex);
		obs_t obs = obss[epoIndex];
		/* 参数个数 xyzt 四个+1个ZTD+2个ZTD梯度+卫星个数*模糊度 */

		if (time.time != 0)rtk.tt = timediff(rtk.sol.time, time);
		udstate_ppp(obs);
  
		int NV = obs.nmax * 2;
		int NX = rtk.nx;
		double dr[3] = { 0 };
		tidedisp(gpst2utc(obs.data[0].time), rtk.sol.rr, odisp, dr);

		CMatrix x1, P1, dx;
		CMatrix H(NV, NX), R(NV, NV), v(NV, 1), I(NX, NX);
		CMatrix Q, K;
		R.Unit();

		x1 = rtk.x;
		P1 = rtk.P;

		/* ppp residual */
		FormKFDesConsMatEachEpoch(0, epoIndex, obs, SatPosClk, x1, H, R, v, dr);

		/* Kalman Filter */
		int* iv = new int[NV], * ix = new int[NX];
		int i, j;
		int nx = 0, nv = 0;
		for (i = nv = 0; i < NV; i++)if (H(i, 0) != 0)iv[nv++] = i;
		for (i = nx = 0; i < NX; i++)if (x1(i, 0) != 0 && P1(i, i) > 0)ix[nx++] = i;
		//xc:para counts(!=0),vc: residual counts(!=0);
		CMatrix H_(nv, nx), R_(nv, nv), v_(nv, 1), x_(nx, 1), P_(nx, nx), I_(nx, nx);
		I_.Unit();

		for (i = 0; i < nx; i++)//col
		{
			x_(i, 0) = x1(ix[i], 0);
			for (j = 0; j < nx; j++)P_(j, i) = P1(ix[j], ix[i]);
			for (j = 0; j < nv; j++)//row
			{
				H_(j, i) = H(iv[j], ix[i]);
				R_(j, j) = R(iv[j], iv[j]);
				v_(j, 0) = v(iv[j], 0);
			}
		}


		CMatrix HT_ = ~H_;
		Q = H_ * P_ * HT_ + R_;
		CMatrix QInv = Q.Inv();
		K = P_ * HT_ * QInv;
	 	CMatrix x1_ = x_ + K * v_;
		CMatrix dx_ =  K * v_;
		CMatrix P1_ = (I_ - K * H_) * P_;

		//PrintMatrix(P1_, testMat);
			
		for (i = 0; i < nx; i++)//col
		{
			x1(ix[i], 0) = x1_(i, 0);
			for (j = 0; j < nx; j++)P1(ix[i], ix[j]) = P1_(i, j);
		}
			
		rtk.x = x1;
		rtk.P = P1;

		for (int i = 0; i < 3; i++)
		{
			sol.data[epoIndex].rr[i] = rtk.x(i, 0);
			sol.data[epoIndex].qr[i] = rtk.P(i, i);
		}
		sol.data[epoIndex].dtr = rtk.x(3, 0) / CLIGHT;
		sol.data[epoIndex].qr[3] = rtk.P(0, 1);
		sol.data[epoIndex].qr[4] = rtk.P(0, 2);
		sol.data[epoIndex].qr[5] = rtk.P(1, 2);
	}
	double a = 1;
}



void CBdPPP::udpos_ppp()
{
	/* init first epoch */
	if (rtk.x(0, 0) == 0 && rtk.x(1, 0) == 0 && rtk.x(2, 0) == 0)
	{
		for (int i = 0; i < 3; i++)
		{
			rtk.x(i, 0) = rtk.sol.rr[i];
			rtk.P(i, i) = VAR_POS;
		}
	}
	/* update else epoch */
	for (int i = 0; i < 3; i++) 
	{
		rtk.P(i, i) += SQR(prn[5]) * fabs(rtk.tt);
	}
}

void CBdPPP::udclk_ppp()
{
	rtk.x(3, 0) = rtk.sol.dtr * CLIGHT;
	rtk.P(3, 3) = VAR_CLK;
	for (int i = 0; i < rtk.nx; i++)if (i != 3)rtk.P(3, i) = rtk.P(i, 3) = 0;
}

void CBdPPP::udtrop_ppp()
{
	double blh[3], azel[] = { 0.0,PI / 2.0 }, ztd, var;

	/* init first epoch */
	if (rtk.x(4,0) == 0.0) 
	{
		ecef2blh(rtk.sol.rr, blh);
		ztd = sbstropcorr(rtk.sol.time, blh, azel, &var);
		rtk.x(4, 0) = ztd;
		rtk.P(4, 4) = var;
		for (int i = 5; i < 7; i++) 
		{
			rtk.x(i, 0) = 1e-6;
			rtk.P(i, i) = VAR_GRA;
		}
	}
	/* update else epoch */
	else 
	{
		rtk.P(4, 4) += SQR(prn[2]) * fabs(rtk.tt);
		for (int i = 5; i < 7; i++) rtk.P(i, i) += SQR(prn[2]*0.1) * fabs(rtk.tt);
	}
}

void CBdPPP::udbias_ppp(obs_t obs)
{
	int nobs = obs.nmax;
	double L[NFREQ], P[NFREQ], Lc, Pc, *bias, offset = 0.0, blh[3] = { 0 };
	double freq1, freq2, ion, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	int i, j, k, f, sat, *slip, clk_jump = 0;
	bias = new double[nobs];
	slip = new int[nobs];

	/* init slip state of all sats */
	for (int i = 0; i < NSATCMP; i++)
		for (int j = 0; j < NFREQ; j++)
			rtk.ssat[i].slip[j] = 0;

	/* detect phase slip */
	detslp_gf(obs);
	detslp_mw(obs);

	ecef2blh(rtk.sol.rr, blh);

	//////////////k记录无周跳观测值个数
	for (i = k = 0; i < nobs ; i++)
	{

		j = obs.data[i].prn - 1 + 7;
		if (i == 10)
		{
			double a = 1;
		}
		corr_meas(obs.data[i], rtk.ssat[obs.data[i].prn - 1].azel, dantr, dants, 0.0, L, P, &Lc, &Pc);
		bias[i] = 0.0;
		bias[i] = Lc - Pc;
		slip[i] = rtk.ssat[obs.data[i].prn - 1].slip[0] || rtk.ssat[obs.data[i].prn - 1].slip[1];
		if (rtk.x(j,0) == 0.0 || slip[i] || bias[i] == 0.0) continue;
		offset += bias[i] - rtk.x(j, 0);
		k++;
	}
	/* correct phase-code jump to ensure phase-code coherency */
	if (k >= 2 && fabs(offset / k) > 0.0005 * CLIGHT) 
	{
		for (i = 0; i < nobs; i++) 
		{
			j = obs.data[i].prn - 1 + 7;
			if (rtk.x(j, 0) != 0.0) rtk.x(j, 0) += offset / k;
		}
	}
	for (i = 0; i < nobs; i++) 
	{
		j = obs.data[i].prn - 1 + 7;
		rtk.P(j,j) += SQR(prn[0]) * fabs(rtk.tt);
		if (bias[i] == 0.0 || (rtk.x(j, 0) != 0.0 && !slip[i])) continue;

		/* reinitialize phase-bias if detecting cycle slip */
		rtk.x(j, 0) = bias[i];
		rtk.P(j, j) = VAR_BIAS;
		for (int i = 0; i < rtk.nx; i++)if (i != j)rtk.P(j, i) = rtk.P(i, j) = 0;
	}
	delete[]slip;
	slip = NULL;
}

double CBdPPP::sbstropcorr(gtime_t time, const double* blh, const double* azel, double* var)
{
	const double k1 = 77.604, k2 = 382000.0, rd = 287.054, gm = 9.784, g = 9.80665;
	static double blh_[3] = { 0 }, zh = 0.0, zw = 0.0;
	int i;
	double c, met[10], sinel = sin(azel[1]), h = blh[2], m;

	if (blh[2] < -100.0 || 10000.0 < blh[2] || azel[1] <= 0) {
		*var = 0.0;
		return 0.0;
	}
	if (zh == 0.0 || fabs(blh[0] - blh_[0]) > 1E-7 || fabs(blh[1] - blh_[1]) > 1E-7 ||
		fabs(blh[2] - blh_[2]) > 1.0) {
		getmet(blh[0] * 180.0 / PI, met);
		c = cos(2.0 * PI * (time2doy(time) - (blh[0] >= 0.0 ? 28.0 : 211.0)) / 365.25);
		for (i = 0; i < 5; i++) met[i] -= met[i + 5] * c;
		zh = 1E-6 * k1 * rd * met[0] / gm;
		zw = 1E-6 * k2 * rd / (gm * (met[4] + 1.0) - met[3] * rd) * met[2] / met[1];
		zh *= pow(1.0 - met[3] * h / met[1], g / (rd * met[3]));
		zw *= pow(1.0 - met[3] * h / met[1], (met[4] + 1.0) * g / (rd * met[3]) - 1.0);
		for (i = 0; i < 3; i++) blh_[i] = blh[i];
	}
	m = 1.001 / sqrt(0.002001 + sinel * sinel);
	*var = 0.12 * 0.12 * m * m;
	return (zh + zw) * m;
}

double CBdPPP::time2doy(gtime_t t)
{
	double ep[6];

	time2epoch(t, ep);
	ep[1] = ep[2] = 1.0; ep[3] = ep[4] = ep[5] = 0.0;
	return timediff(t, epoch2time(ep)) / 86400.0 + 1.0;
}

gtime_t CBdPPP::epoch2time(const double* ep)
{
	const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
	gtime_t time = { 0 };
	int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];

	if (year < 1970 || 2099 < year || mon < 1 || 12 < mon) return time;

	/* leap year if year%4==0 in 1901-2099 */
	days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
	sec = (int)floor(ep[5]);
	time.time = (time_t)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
	time.sec = ep[5] - sec;
	return time;
}

gtime_t CBdPPP::utc2gpst(gtime_t t)
{
	int i;

	for (i = 0; leaps[i][0] > 0; i++) {
		if (timediff(t, epoch2time(leaps[i])) >= 0.0) return timeadd(t, -leaps[i][6]);
	}
	return t;
}

gtime_t CBdPPP::gpst2utc(gtime_t t)
{
	gtime_t tu;
	int i;

	for (i = 0; leaps[i][0] > 0; i++) {
		tu = timeadd(t, leaps[i][6]);
		if (timediff(tu, epoch2time(leaps[i])) >= 0.0) return tu;
	}
	return t;
}

double CBdPPP::utc2gmst(gtime_t t, double ut1_utc)
{
	const double ep2000[] = { 2000,1,1,12,0,0 };
	gtime_t tut, tut0;
	double ut, t1, t2, t3, gmst0, gmst;

	tut = timeadd(t, ut1_utc);
	ut = time2sec(tut, &tut0);
	t1 = timediff(tut0, epoch2time(ep2000)) / 86400.0 / 36525.0;
	t2 = t1 * t1; t3 = t2 * t1;
	gmst0 = 24110.54841 + 8640184.812866 * t1 + 0.093104 * t2 - 6.2E-6 * t3;
	gmst = gmst0 + 1.002737909350795 * ut;

	return fmod(gmst, 86400.0) * PI / 43200.0; /* 0 <= gmst <= 2*PI */
}

double CBdPPP::time2sec(gtime_t time, gtime_t* day)
{
	double ep[6], sec;
	time2epoch(time, ep);
	sec = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
	ep[3] = ep[4] = ep[5] = 0.0;
	*day = epoch2time(ep);
	return sec;
}

void CBdPPP::getmet(double lat, double* met)
{
	static const double metprm[][10] = { /* lat=15,30,45,60,75 */
		{1013.25,299.65,26.31,6.30E-3,2.77,  0.00, 0.00,0.00,0.00E-3,0.00},
		{1017.25,294.15,21.79,6.05E-3,3.15, -3.75, 7.00,8.85,0.25E-3,0.33},
		{1015.75,283.15,11.66,5.58E-3,2.57, -2.25,11.00,7.24,0.32E-3,0.46},
		{1011.75,272.15, 6.78,5.39E-3,1.81, -1.75,15.00,5.36,0.81E-3,0.74},
		{1013.00,263.65, 4.11,4.53E-3,1.55, -0.50,14.50,3.39,0.62E-3,0.30}
	};
	int i, j;
	double a;
	lat = fabs(lat);
	if (lat <= 15.0) for (i = 0; i < 10; i++) met[i] = metprm[0][i];
	else if (lat >= 75.0) for (i = 0; i < 10; i++) met[i] = metprm[4][i];
	else {
		j = (int)(lat / 15.0); a = (lat - j * 15.0) / 15.0;
		for (i = 0; i < 10; i++) met[i] = (1.0 - a) * metprm[j - 1][i] + a * metprm[j][i];
	}
}

void CBdPPP::detslp_gf(obs_t obs)
{
	double g0, g1;
	int i, j;
	int n = obs.nmax;
	for (i = 0; i < n ; i++) 
	{
		if ((g1 = gfmeas(obs.data[i])) == 0.0) continue;

		g0 = rtk.ssat[obs.data[i].prn - 1].gf[0];
		rtk.ssat[obs.data[i].prn - 1].gf[0] = g1;

		if (g0 != 0.0 && fabs(g1 - g0) > thresslip) 
		{
			for (j = 0; j < NFREQ; j++) rtk.ssat[obs.data[i].prn - 1].slip[j] =1;
		}
	}
}

void CBdPPP::detslp_mw(obs_t obs)
{
	double w0, w1;
	int i, j;
	int n = obs.nmax;

	for (i = 0; i < n ; i++) {
		if ((w1 = mwmeas(obs.data[i])) == 0.0) continue;

		w0 = rtk.ssat[obs.data[i].prn - 1].mw[0];
		rtk.ssat[obs.data[i].prn - 1].mw[0] = w1;

		if (w0 != 0.0 && fabs(w1 - w0) > THRES_MW_JUMP) 
		{
			for (j = 0; j < NFREQ; j++)rtk.ssat[obs.data[i].prn - 1].slip[j] = 1;
		}
	}
}

double CBdPPP::gfmeas(obsd_t obs)
{
	if (obs.L[0] == 0.0 || obs.L[1] == 0.0) return 0.0;
	return (obs.L[0] / FREQ1_CMP - obs.L[1] / FREQ3_CMP) * CLIGHT;
}

double CBdPPP::mwmeas(obsd_t obs)
{
	if ( obs.L[0] == 0.0 || obs.L[1] == 0.0 ||
		obs.P[0] == 0.0 || obs.P[1] == 0.0) return 0.0;
	return (obs.L[0] - obs.L[1]) * CLIGHT / (FREQ1_CMP - FREQ3_CMP) -
		(FREQ1_CMP * obs.P[0] + FREQ3_CMP * obs.P[1]) / (FREQ1_CMP + FREQ3_CMP);
}

 void CBdPPP::corr_meas(obsd_t obs, const double* azel, const double* dantr, const double* dants, double phw, double* L, double* P, double* Lc, double* Pc)
{
	double lam[2];
	double C1, C2;
	int i, sys;
	int prn; /////////PANLIN
	int flag_row;//////////PANLIN
	double el_deg;//////////PANLIN
	dcb_t cbias;

	/* search for dcb correction data of this sat */
	cbias.bias = 0;
	for (int ind = 0; ind < MAXSAT; ind++)
	{
		if (navs.cbias[ind].sat == obs.sat)
		{
			cbias = navs.cbias[ind];
			break;
		}
	}

	lam[0] = CLIGHT / FREQ1_CMP;
	lam[1] = CLIGHT / FREQ3_CMP;

	for (i = 0; i < NFREQ; i++)
	{
		L[i] = P[i] = 0.0;
		if (lam[i] == 0.0 || obs.L[i] == 0.0 || obs.P[i] == 0.0) continue;

		/* antenna phase center and phase windup correction */
		L[i] = obs.L[i] * lam[i] - dants[i] - dantr[i] - phw * lam[i];
		P[i] = obs.P[i] - dants[i] - dantr[i];

		/* C2I-C6I dcb correction */
		/*
		if (i == 1)
		{
			P[i] += cbias.bias;
		}
		*/
#if 0
		L[i] -= 0.25 * lam[i]; *1 / 4 cycle - shift *
#endif
	}
		

	/* iono-free LC */
	*Lc = *Pc = 0.0;
	if (lam[0] == 0.0 || lam[1] == 0.0) return;

	C1 = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
	C2 = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));

#if 1  ///////////////////////////////////////////////////////////PANLIN   change 0 to 1 for uncombined PPP
	/* C2I-C6I dcb correction (C2I->Pc,C6I->Pc) */
	{////////////////////////////add    |SYS_GAL|SYS_CMP
		/*
		if (P[0] != 0.0) P[0] -= C2 * cbias.bias;
		if (P[1] != 0.0) P[1] += C1 * cbias.bias;           ////////////////////////    if (P[1]!=0.0) P[1]+=C1*nav->cbias[obs->sat-1][0];
		*/
	}
#endif
	///////////////////////////SICB correction   PANLIN
	/////////////+SICB???  -SICB???
	el_deg = azel[1] * 180.0 / PI;
	if ((el_deg > 0.0) && (el_deg < 90.0))
	{
		flag_row = floor(el_deg / 10.0);

		/////IGSO
		if ((obs.sat == _T("C06")) || (obs.sat == _T("C07")) || (obs.sat == _T("C08")) || (obs.sat == _T("C09")) || (obs.sat == _T("C10")) || (obs.sat == _T("C13")) || (obs.sat == _T("C16")))
		{
			if (P[0] != 0.0) P[0] = P[0] + (BDS_IGSO[0][flag_row] + (BDS_IGSO[0][flag_row + 1] - BDS_IGSO[0][flag_row]) * (el_deg - flag_row * 10.0) / 10.0);
			if (P[i] != 0.0) P[i] = P[i] + (BDS_IGSO[i][flag_row] + (BDS_IGSO[i][flag_row + 1] - BDS_IGSO[i][flag_row]) * (el_deg - flag_row * 10.0) / 10.0);
		}

		/////MEO
		if ((obs.sat == _T("C11")) || (obs.sat == _T("C12")) || (obs.sat == _T("C14")))
		{
			if (P[0] != 0.0) P[0] = P[0] + (BDS_MEO[0][flag_row] + (BDS_MEO[0][flag_row + 1] - BDS_MEO[0][flag_row]) * (el_deg - flag_row * 10.0) / 10.0);
			if (P[i] != 0.0) P[i] = P[i] + (BDS_MEO[i][flag_row] + (BDS_MEO[i][flag_row + 1] - BDS_MEO[i][flag_row]) * (el_deg - flag_row * 10.0) / 10.0);
		}
	}
	///////////////////////////SICB correction

	if (L[0] != 0.0 && L[1] != 0.0) *Lc = C1 * L[0] + C2 * L[1];
	if (P[0] != 0.0 && P[1] != 0.0) *Pc = C1 * P[0] + C2 * P[1];
}

double CBdPPP::interppol(const double* x, double* y, int n)
 {
	 int i, j;

	 for (j = 1; j < n; j++) 
	 {
		 for (i = 0; i < n - j; i++) 
		 {
			 y[i] = (x[i + j] * y[i] - x[i] * y[i + 1]) / (x[i + j] - x[i]);
		 }
	 }
	 return y[0];
 }

int CBdPPP::pephpos(double* rs, gtime_t time, int epoIndex, int satIndex)
{
	int i, j, k, index, flag = 0;
	double t[NMAX + 1], p[3][NMAX + 1], c[2], * pos, std = 0.0, s[3], sinl, cosl;

	/* binary search */
	for (i = 0, j = navs.np - 1; i < j;)
	{
		k = (i + j) / 2;
		if (timediff(navs.peph[k].time, time) < 0.0) i = k + 1; else j = k;
	}
	index = i <= 0 ? 0 : i - 1;

	/* polynomial interpolation for orbit */
	i = index - (NMAX + 1) / 2;
	if (i < 0) i = 0; else if (i + NMAX >= navs.np) i = navs.np - NMAX - 1;
	/* time of sat moving from peph time */
	for (j = 0; j <= NMAX; j++)
	{
		t[j] = timediff(navs.peph[i + j].time, time);
	}
	/* correciton for earh rotation ver.2.4.0 */
	for (j = 0; j <= NMAX; j++)
	{
		for (int satInd = 0; satInd < navs.peph[i + j].SatCount; satInd++)
		{
			if (obss[epoIndex].data[satIndex].sat == navs.peph[i + j].satstate[satInd].sat)
			{
				pos = navs.peph[i + j].satstate[satInd].pos;
				sinl = sin(OMGe * t[j]);
				cosl = cos(OMGe * t[j]);
				p[0][j] = cosl * pos[0] - sinl * pos[1];
				p[1][j] = sinl * pos[0] + cosl * pos[1];
				p[2][j] = pos[2];
				flag = 1;//if found the peph of this sat;
				break;
			}
		}
	}
	if (!flag)return 0;
	/* interpolation */
	for (i = 0; i < 3; i++)
	{
		rs[i] = interppol(t, p[i], NMAX + 1);
	}
	return 1;
}

int CBdPPP::pephclk(double& dts, gtime_t time, int epoIndex, int satIndex)
{
	double t[2], c[2], std;
	int i, j, k, index, flag = 0;
	/* binary search */
	for (i = 0, j = navs.nc - 1; i < j;) {
		k = (i + j) / 2;
		if (timediff(navs.pclks[k].time, time) < 0.0) i = k + 1; else j = k;
	}
	index = i <= 0 ? 0 : i - 1;

	/* linear interpolation for clock */
	t[0] = timediff(time, navs.pclks[index].time);
	t[1] = timediff(time, navs.pclks[index + 1].time);
	for (int satInd = 0; satInd < navs.pclks[index].nmax; satInd++)
	{
		if (obss[epoIndex].data[satIndex].sat == navs.pclks[index].pclk[satInd].sat)
		{
			c[0] = navs.pclks[index].pclk[satInd].clk;
			c[1] = navs.pclks[index + 1].pclk[satInd].clk;
			flag = 1;
			break;
		}
	}
	if (!flag)return 0;
	if (t[0] <= 0.0)
	{
		if ((dts = c[0]) == 0.0) return 0;
		//std = nav->pclk[index].std[sat - 1][0] * CLIGHT - EXTERR_CLK * t[0];
	}
	else if (t[1] >= 0.0)
	{
		if ((dts = c[1]) == 0.0) return 0;
		//std = nav->pclk[index + 1].std[sat - 1][0] * CLIGHT + EXTERR_CLK * t[1];
	}
	else if (c[0] != 0.0 && c[1] != 0.0)
	{
		dts = (c[1] * t[0] - c[0] * t[1]) / (t[0] - t[1]);
		i = t[0] < -t[1] ? 0 : 1;
		//std = nav->pclk[index + i].std[sat - 1][0] * CLIGHT + EXTERR_CLK * fabs(t[i]);
	}
	return 1;
}

void CBdPPP::satantoff(gtime_t time, const double* rs, CString sat, double* dant)
{
	pcv_t pcv;
	double rsun[3], gmst, erpv[5] = { 0 };
	long double ex[3], ey[3], ez[3], es[3], r[3], freq[2], norm;
	double C1, C2, dant1, dant2;
	int i;

	for (i = 0; i < MAXSAT; i++) 
	if (navs.pcvs[i].sat == sat && (timediff(time, navs.pcvs[i].ts) > 0 && timediff(time, navs.pcvs[i].te) < 0))
	{
		pcv = navs.pcvs[i];
		break;
	}
	dant[0] = dant[1] = dant[2] = 0.0;

	/* sun position in ecef */
	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, &gmst);

	/* unit vectors of satellite fixed coordinates */

	for (i = 0; i < 3; i++) r[i] = -rs[i];
	norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
	for (i = 0; i < 3; i++)ez[i] = r[i] / norm;

	for (i = 0; i < 3; i++) r[i] = rsun[i] - rs[i];
	norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
	for (i = 0; i < 3; i++)es[i] = r[i] / norm;

	{
		r[0] = ez[1] * es[2] - ez[2] * es[1];
		r[1] = ez[2] * es[0] - ez[0] * es[2];
		r[2] = ez[0] * es[1] - ez[1] * es[0];
	}

	norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
	for (i = 0; i < 3; i++)ey[i] = r[i] / norm;

	{
		ex[0] = ey[1] * ez[2] - ey[2] * ez[1];
		ex[1] = ey[2] * ez[0] - ey[0] * ez[2];
		ex[2] = ey[0] * ez[1] - ey[1] * ez[0];
	}

	/* iono-free LC coefficients */
	freq[0] = FREQ1_CMP;
	freq[1] = FREQ3_CMP;

	C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
	C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

	/* iono-free LC */
	for (i = 0; i < 3; i++) {
		dant1 = pcv.off[0][0] * ex[i] + pcv.off[0][1] * ey[i] + pcv.off[0][2] * ez[i];
		dant2 = pcv.off[1][0] * ex[i] + pcv.off[1][1] * ey[i] + pcv.off[1][2] * ez[i];
		dant[i] = C1 * dant1 + C2 * dant2;
	}
}

void CBdPPP::sunmoonpos(gtime_t tutc, const double* erpv, double* rsun, double* rmoon, double* gmst)
{
	gtime_t tut;
	double rs[3], rm[3], U[9], gmst_;

	tut = timeadd(tutc, erpv[2]); /* utc -> ut1 */

	/* sun and moon position in eci */
	sunmoonpos_eci(tut, rsun ? rs : NULL, rmoon ? rm : NULL);

	/* eci to ecef transformation matrix */
	eci2ecef(tutc, erpv, U, &gmst_);

	CMatrix Umat(3,3), rsmat(3, 1), rsunmat(3, 1), rmmat(3, 1), rmoonmat(3, 1);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) Umat(i, j) = U[i + j * 3];
		if (rsun)rsmat(i, 0) = rs[i];
		if (rmoon)rmmat(i, 0) = rm[i];
	}

	/* sun and moon postion in ecef */
	if (rsun) rsunmat = Umat * rsmat;
	if (rmoon) rmoonmat = Umat * rmmat;
	if (gmst) *gmst = gmst_;
	for (int i = 0; i < 3; i++)
	{
		if (rsun)rsun[i] = rsunmat(i, 0);
		if (rmoon)rmoon[i] = rmoonmat(i, 0);
	}
}

void CBdPPP::sunmoonpos_eci(const gtime_t tut, double* rsun, double* rmoon)
{
	CString strEp2000[] = { _T("2000"),_T("1"),_T("1"),_T("12"),_T("0"),_T("0") };
	double t, f[5], eps, Ms, ls, rs, lm, pm, rm, sine, cose, sinp, cosp, sinl, cosl;

	gtime_t ep2000;
	str2time(strEp2000, ep2000);
	t = timediff(tut, ep2000) / 86400.0 / 36525.0;

	/* astronomical arguments */
	ast_args(t, f);

	/* obliquity of the ecliptic */
	eps = 23.439291 - 0.0130042 * t;
	sine = sin(eps * PI / 180.0); cose = cos(eps * PI / 180.0);

	/* sun position in eci */
	if (rsun) {
		Ms = 357.5277233 + 35999.05034 * t;
		ls = 280.460 + 36000.770 * t + 1.914666471 * sin(Ms * PI / 180.0) + 0.019994643 * sin(2.0 * Ms * PI / 180.0);
		rs = AU * (1.000140612 - 0.016708617 * cos(Ms * PI / 180.0) - 0.000139589 * cos(2.0 * Ms * PI / 180.0));
		sinl = sin(ls * PI / 180.0); cosl = cos(ls * PI / 180.0);
		rsun[0] = rs * cosl;
		rsun[1] = rs * cose * sinl;
		rsun[2] = rs * sine * sinl;

	}
	/* moon position in eci */
	if (rmoon) {
		lm = 218.32 + 481267.883 * t + 6.29 * sin(f[0]) - 1.27 * sin(f[0] - 2.0 * f[3]) +
			0.66 * sin(2.0 * f[3]) + 0.21 * sin(2.0 * f[0]) - 0.19 * sin(f[1]) - 0.11 * sin(2.0 * f[2]);
		pm = 5.13 * sin(f[2]) + 0.28 * sin(f[0] + f[2]) - 0.28 * sin(f[2] - f[0]) -
			0.17 * sin(f[2] - 2.0 * f[3]);
		rm = RE_WGS84 / sin((0.9508 + 0.0518 * cos(f[0]) + 0.0095 * cos(f[0] - 2.0 * f[3]) +
			0.0078 * cos(2.0 * f[3]) + 0.0028 * cos(2.0 * f[0])) * PI / 180.0);
		sinl = sin(lm * PI / 180.0); cosl = cos(lm * PI / 180.0);
		sinp = sin(pm * PI / 180.0); cosp = cos(pm * PI / 180.0);
		rmoon[0] = rm * cosp * cosl;
		rmoon[1] = rm * (cose * cosp * sinl - sine * sinp);
		rmoon[2] = rm * (sine * cosp * sinl + cose * sinp);

	}
}

void CBdPPP::eci2ecef(gtime_t tutc, const double* erpv, double* U, double* gmst)
{
	const double ep2000[] = { 2000,1,1,12,0,0 };
	static gtime_t tutc_;
	static double U_[9], gmst_;
	gtime_t tgps;
	double eps, ze, th, z, t, t2, t3, dpsi, deps, gast, f[5];
	double R1[9], R2[9], R3[9], R[9], W[9], N[9], P[9], NP[9];
	int i;

	if (fabs(timediff(tutc, tutc_)) < 0.01) { /* read cache */
		for (i = 0; i < 9; i++) U[i] = U_[i];
		if (gmst) *gmst = gmst_;
		return;
	}
	tutc_ = tutc;

	/* terrestrial time */
	tgps = utc2gpst(tutc_);
	t = (timediff(tgps, epoch2time(ep2000)) + 19.0 + 32.184) / 86400.0 / 36525.0;
	t2 = t * t; t3 = t2 * t;

	/* astronomical arguments */
	ast_args(t, f);

	/* iau 1976 precession */
	ze = (2306.2181 * t + 0.30188 * t2 + 0.017998 * t3) * PI / 180.0 / 3600.0;
	th = (2004.3109 * t - 0.42665 * t2 - 0.041833 * t3) * PI / 180.0 / 3600.0;
	z = (2306.2181 * t + 1.09468 * t2 + 0.018203 * t3) * PI / 180.0 / 3600.0;
	eps = (84381.448 - 46.8150 * t - 0.00059 * t2 + 0.001813 * t3) * PI / 180.0 / 3600.0;
	Rz(-z, R1); Ry(th, R2); Rz(-ze, R3);

	matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
	matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

	/* iau 1980 nutation */
	nut_iau1980(t, f, &dpsi, &deps);
	Rx(-eps - deps, R1); Rz(-dpsi, R2); Rx(eps, R3);
	matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
	matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

	/* greenwich aparent sidereal time (rad) */
	gmst_ = utc2gmst(tutc_, erpv[2]);
	gast = gmst_ + dpsi * cos(eps);
	gast += (0.00264 * sin(f[4]) + 0.000063 * sin(2.0 * f[4])) * PI / 180.0 / 3600;

	/* eci to ecef transformation matrix */
	Ry(-erpv[0], R1); Rx(-erpv[1], R2); Rz(gast, R3);
	matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, W);
	matmul("NN", 3, 3, 3, 1.0, W, R3, 0.0, R); /* W=Ry(-xp)*Rx(-yp) */
	matmul("NN", 3, 3, 3, 1.0, N, P, 0.0, NP);
	matmul("NN", 3, 3, 3, 1.0, R, NP, 0.0, U_); /* U=W*Rz(gast)*N*P */

	for (i = 0; i < 9; i++) U[i] = U_[i];
	if (gmst) *gmst = gmst_;

}

void CBdPPP::ast_args(double t, double* f)
{
	static const double fc[][5] = { /* coefficients for iau 1980 nutation */
		{ 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470},
		{ 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149},
		{  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417},
		{ 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169},
		{ 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939}
	};
	double tt[4];
	int i, j;

	for (tt[0] = t, i = 1; i < 4; i++) tt[i] = tt[i - 1] * t;
	for (i = 0; i < 5; i++) {
		f[i] = fc[i][0] * 3600.0;
		for (j = 0; j < 4; j++) f[i] += fc[i][j + 1] * tt[j];
		f[i] = fmod(f[i] * PI/180.0/3600.0, 2.0 * PI);
	}
}

void CBdPPP::nut_iau1980(double t, const double* f, double* dpsi, double* deps)
{
	static const double nut[106][10] = {
		{   0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9},
		{   0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1},
		{   0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5},
		{   0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5},
		{   0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1},
		{   1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0},
		{   0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6},
		{   0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0},
		{   1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1},
		{   0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3},
		{  -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0},
		{   0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0},
		{  -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0},
		{   1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0},
		{   0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0},
		{  -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0},
		{  -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0},
		{   1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0},
		{  -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0},
		{  -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0},
		{   0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0},
		{   2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0},
		{   2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0},
		{   1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0},
		{   0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0},
		{   0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0},
		{  -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0},
		{   0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0},
		{   0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0},
		{  -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0},
		{   0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0},
		{   1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0},
		{   0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0},
		{   2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0},
		{  -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0},
		{   1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0},
		{   0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0},
		{   0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0},
		{   1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0},
		{   0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0},
		{  -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0},
		{   0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0},
		{   2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0},
		{   1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0},
		{   1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0},
		{   0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0},
		{   0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0},
		{   2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0},
		{   1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0},
		{   1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0},
		{   0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0},
		{   0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0},
		{   1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0},
		{   2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0},
		{   0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0},
		{   1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0},
		{   1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0},
		{  -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0},
		{   0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0},
		{   1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0},
		{   3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0},
		{  -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0},
		{   1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0},
		{  -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0},
		{   1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0},
		{  -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0},
		{   0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0},
		{  -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0},
		{   2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0},
		{   3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0},
		{   1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0},
		{   0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0},
		{   1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0},
		{   1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0},
		{   1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0},
		{   0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0},
		{   0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0},
		{   0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0},
		{   1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0},
		{   1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0},
		{   1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0},
		{   1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0},
		{   2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0},
		{   0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0},
		{   0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0},
		{  -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0},
		{   2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0},
		{   0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0},
		{   0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0},
		{   0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0},
		{   0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0},
		{   1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0},
		{   3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0},
		{  -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0},
		{  -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0},
		{   0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0},
		{   0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0},
		{  -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0},
		{   2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0},
		{   2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0},
		{   2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0},
		{   2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0},
		{   1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0},
		{  -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0},
		{  -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0},
		{   0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0}
	};
	double ang;
	int i, j;

	*dpsi = *deps = 0.0;

	for (i = 0; i < 106; i++) {
		ang = 0.0;
		for (j = 0; j < 5; j++) ang += nut[i][j] * f[j];
		*dpsi += (nut[i][6] + nut[i][7] * t) * sin(ang);
		*deps += (nut[i][8] + nut[i][9] * t) * cos(ang);
	}
	*dpsi *= 1E-4 * PI / 180.0 / 3600; /* 0.1 mas -> rad */
	*deps *= 1E-4 * PI / 180.0 / 3600;
}

void CBdPPP::tide_solid(const double* rsun, const double* rmoon,const double* pos, const double* E, double gmst, double* dr)
{
	double dr1[3], dr2[3], eu[3], du, dn, sinl, sin2l;

	/* step1: time domain */
	eu[0] = E[2]; eu[1] = E[5]; eu[2] = E[8];
	tide_pl(eu, rsun, GMS, pos, dr1);
	tide_pl(eu, rmoon, GMM, pos, dr2);

	/* step2: frequency domain, only K1 radial */
	sin2l = sin(2.0 * pos[0]);
	du = -0.012 * sin2l * sin(gmst + pos[1]);

	dr[0] = dr1[0] + dr2[0] + du * E[2];
	dr[1] = dr1[1] + dr2[1] + du * E[5];
	dr[2] = dr1[2] + dr2[2] + du * E[8];

}

void CBdPPP::tide_oload(gtime_t tut, const double* odisp, double* denu)
{
	const double args[][5] = {
		{1.40519E-4, 2.0,-2.0, 0.0, 0.00},  /* M2 */
		{1.45444E-4, 0.0, 0.0, 0.0, 0.00},  /* S2 */
		{1.37880E-4, 2.0,-3.0, 1.0, 0.00},  /* N2 */
		{1.45842E-4, 2.0, 0.0, 0.0, 0.00},  /* K2 */
		{0.72921E-4, 1.0, 0.0, 0.0, 0.25},  /* K1 */
		{0.67598E-4, 1.0,-2.0, 0.0,-0.25},  /* O1 */
		{0.72523E-4,-1.0, 0.0, 0.0,-0.25},  /* P1 */
		{0.64959E-4, 1.0,-3.0, 1.0,-0.25},  /* Q1 */
		{0.53234E-5, 0.0, 2.0, 0.0, 0.00},  /* Mf */
		{0.26392E-5, 0.0, 1.0,-1.0, 0.00},  /* Mm */
		{0.03982E-5, 2.0, 0.0, 0.0, 0.00}   /* Ssa */
	};
	const double ep1975[] = { 1975,1,1,0,0,0 };
	double ep[6], fday, days, t, t2, t3, a[5], ang, dp[3] = { 0 };
	int i, j;

	/* angular argument: see subroutine arg.f for reference [1] */
	time2epoch(tut, ep);
	fday = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
	ep[3] = ep[4] = ep[5] = 0.0;
	days = timediff(epoch2time(ep), epoch2time(ep1975)) / 86400.0 + 1.0;
	t = (27392.500528 + 1.000000035 * days) / 36525.0;
	t2 = t * t; t3 = t2 * t;

	a[0] = fday;
	a[1] = (279.69668 + 36000.768930485 * t + 3.03E-4 * t2) * PI / 180.0;/* H0 */
	a[2] = (270.434358 + 481267.88314137 * t - 0.001133 * t2 + 1.9E-6 * t3) * PI / 180.0; /* S0 */
	a[3] = (334.329653 + 4069.0340329577 * t - 0.010325 * t2 - 1.2E-5 * t3) * PI / 180.0; /* P0 */
	a[4] = 2.0 * PI;

	/* displacements by 11 constituents */
	for (i = 0; i < 11; i++) {
		ang = 0.0;
		for (j = 0; j < 5; j++) ang += a[j] * args[i][j];
		for (j = 0; j < 3; j++) dp[j] += odisp[j + i * 6] * cos(ang - odisp[j + 3 + i * 6] * PI / 180.0);
	}
	denu[0] = -dp[1];
	denu[1] = -dp[2];
	denu[2] = dp[0];

}

void CBdPPP::tide_pole(gtime_t tut, const double* pos, const double* erpv, double* denu)
{
	double xp_bar, yp_bar, m1, m2, cosl, sinl;

	/* iers mean pole (mas) */
	iers_mean_pole(tut, &xp_bar, &yp_bar);

	/* ref [7] eq.7.24 */
	m1 = erpv[0] / (PI / 180.0 / 3600.0) - xp_bar * 1E-3; /* (as) */
	m2 = -erpv[1] / (PI / 180.0 / 3600.0) + yp_bar * 1E-3;

	/* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
	cosl = cos(pos[1]);
	sinl = sin(pos[1]);
	denu[0] = 9E-3 * sin(pos[0]) * (m1 * sinl - m2 * cosl); /* de= Slambda (m) */
	denu[1] = -9E-3 * cos(2.0 * pos[0]) * (m1 * cosl + m2 * sinl); /* dn=-Stheta  (m) */
	denu[2] = -33E-3 * sin(2.0 * pos[0]) * (m1 * cosl + m2 * sinl); /* du= Sr      (m) */

}

void CBdPPP::tide_pl(const double* eu, const double* rp, double GMp, const double* pos, double* dr)
{
	const double H3 = 0.292, L3 = 0.015;
	double r, ep[3], latp, lonp, p, K2, K3, a, H2, L2, dp, du, cosp, sinl, cosl;
	int i;

    r = sqrt(pow(rp[0], 2) + pow(rp[1], 2) + pow(rp[2], 2));
	for (i = 0; i < 3; i++) ep[i] = rp[i] / r;

	K2 = GMp / GM * SQR(RE_WGS84) * SQR(RE_WGS84) / (r * r * r);
	K3 = K2 * RE_WGS84 / r;
	latp = asin(ep[2]); lonp = atan2(ep[1], ep[0]);
	cosp = cos(latp); sinl = sin(pos[0]); cosl = cos(pos[0]);

	/* step1 in phase (degree 2) */
	p = (3.0 * sinl * sinl - 1.0) / 2.0;
	H2 = 0.6078 - 0.0006 * p;
	L2 = 0.0847 + 0.0002 * p;
	a = ep[0] * eu[0] + ep[1] * eu[1] + ep[2] * eu[2];
	dp = K2 * 3.0 * L2 * a;
	du = K2 * (H2 * (1.5 * a * a - 0.5) - 3.0 * L2 * a * a);

	/* step1 in phase (degree 3) */
	dp += K3 * L3 * (7.5 * a * a - 1.5);
	du += K3 * (H3 * (2.5 * a * a * a - 1.5 * a) - L3 * (7.5 * a * a - 1.5) * a);

	/* step1 out-of-phase (only radial) */
	du += 3.0 / 4.0 * 0.0025 * K2 * sin(2.0 * latp) * sin(2.0 * pos[0]) * sin(pos[1] - lonp);
	du += 3.0 / 4.0 * 0.0022 * K2 * cosp * cosp * cosl * cosl * sin(2.0 * (pos[1] - lonp));

	dr[0] = dp * ep[0] + du * eu[0];
	dr[1] = dp * ep[1] + du * eu[1];
	dr[2] = dp * ep[2] + du * eu[2];

}

void CBdPPP::iers_mean_pole(gtime_t tut, double* xp_bar, double* yp_bar)
{
	const double ep2000[] = { 2000,1,1,0,0,0 };
	double y, y2, y3;

	y = timediff(tut, epoch2time(ep2000)) / 86400.0 / 365.25;

	if (y < 3653.0 / 365.25) { /* until 2010.0 */
		y2 = y * y; y3 = y2 * y;
		*xp_bar = 55.974 + 1.8243 * y + 0.18413 * y2 + 0.007024 * y3; /* (mas) */
		*yp_bar = 346.346 + 1.7896 * y - 0.10729 * y2 - 0.000908 * y3;
	}
	else { /* after 2010.0 */
		*xp_bar = 23.513 + 7.6141 * y; /* (mas) */
		*yp_bar = 358.891 - 0.6287 * y;
	}
}

void CBdPPP::matmul(const char* tr, int n, int k, int m, double alpha, const double* A, const double* B, double beta, double* C)
{
	double d;
	int i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

	for (i = 0; i < n; i++) for (j = 0; j < k; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x < m; x++) d += A[i + x * n] * B[x + j * m]; break;
		case 2: for (x = 0; x < m; x++) d += A[i + x * n] * B[j + x * k]; break;
		case 3: for (x = 0; x < m; x++) d += A[x + i * m] * B[x + j * m]; break;
		case 4: for (x = 0; x < m; x++) d += A[x + i * m] * B[j + x * k]; break;
		}
		if (beta == 0.0) C[i + j * n] = alpha * d; else C[i + j * n] = alpha * d + beta * C[i + j * n];
	}
}

void CBdPPP::model_trop(gtime_t time, const double* pos, const double* azel, const CMatrix x, double* dtdx, double* dtrp, double* var)
{
	double trp[3] = { 0 }, std[3];
	for (int i = 0; i < 3; i++)trp[i] = x(4 + i, 0);
	*dtrp = trop_model_prec(time, pos, azel, trp, dtdx, var);
}

void CBdPPP::satantpcv(const double* rs, const double* rr, const pcv_t* pcv, double* dant)
{
	double ru[3], rz[3], eu[3], ez[3], nadir, cosa = 0;
	int i;

	for (i = 0; i < 3; i++) {
		ru[i] = rr[i] - rs[i];
		rz[i] = -rs[i];
	}
	for (int i = 0; i < 3; i++)eu[i] = ru[i] / sqrt(pow(ru[0], 2) + pow(ru[1], 2) + pow(ru[2], 2));
	for (int i = 0; i < 3; i++)ez[i] = rz[i] / sqrt(pow(rz[0], 2) + pow(rz[1], 2) + pow(rz[2], 2));

	for (int i = 0; i < 3; i++) cosa += eu[i] * ez[i];
	cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
	nadir = acos(cosa);

	antmodel_s(pcv, nadir, dant);
}

int CBdPPP::model_phw(gtime_t time, CString sat, const double* rs, const double* rr, double* phw)
{
	double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
	double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
	int i;

	/* satellite yaw attitude model */
	if (!sat_yaw(time, sat, rs, exs, eys)) return 0;

	/* unit vector satellite to receiver */
	for (i = 0; i < 3; i++) r[i] = rr[i] - rs[i];
	for (i = 0; i < 3; i++) ek[i] = r[i] / sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));

	/* unit vectors of receiver antenna */
	ecef2blh(rr, pos);
	xyz2enu(pos, E);
	exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
	eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

	/* phase windup effect */
	cross3(ek, eys, eks);
	cross3(ek, eyr, ekr);
	for (i = 0; i < 3; i++) {
		ds[i] = exs[i] - ek[i] * (ek[0] * exs[0] + ek[1] * exs[1] + ek[2] * exs[2]) - eks[i];
		dr[i] = exr[i] - ek[i] * (ek[0] * exr[0] + ek[1] * exr[1] + ek[2] * exr[2]) + ekr[i];
	}
	cosp = (ds[0] * dr[0] + ds[1] * dr[1] + ds[2] * dr[2]) / sqrt(pow(ds[0], 2) + pow(ds[1], 2) + pow(ds[2], 2)) / sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
	if (cosp < -1.0) cosp = -1.0;
	else if (cosp > 1.0) cosp = 1.0;
	ph = acos(cosp) / 2.0 / PI;
	cross3(ds, dr, drs);
	if ((ek[0]* drs[0]+ ek[01] * drs[1]+ ek[2] * drs[2]) < 0.0) ph = -ph;

	*phw = ph + floor(*phw - ph + 0.5); /* in cycle */
	return 1;
}

double CBdPPP::trop_model_prec(gtime_t time, const double* pos, const double* azel, const double* x, double* dtdx, double* var)
{
	const double zazel[] = { 0.0,PI / 2.0 };
	double zhd, m_h, m_w, cotz, grad_n, grad_e;

	/* zenith hydrostatic delay */
	zhd = tropmodel(pos, zazel[1], 0.0);

	/* mapping function */
	m_h = tropmapf(time, pos, azel, &m_w);

	if (azel[1] > 0.0) {

		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
		cotz = 1.0 / tan(azel[1]);
		grad_n = m_w * cotz * cos(azel[0]);
		grad_e = m_w * cotz * sin(azel[0]);
		m_w += grad_n * x[1] + grad_e * x[2];
		dtdx[1] = grad_n * (x[0] - zhd);
		dtdx[2] = grad_e * (x[0] - zhd);
	}
	dtdx[0] = m_w;
	*var = SQR(0.01);
	return m_h * zhd + m_w * (x[0] - zhd);
}

void CBdPPP::antmodel_s(const pcv_t* pcv, double nadir, double* dant)
{
	int i;

	for (i = 0; i < NFREQ; i++) 
	{
		dant[i] = interpvar(nadir * 180.0/PI * 5.0, pcv->var[i]);
	}
}

void CBdPPP::antmodel_r(const pcv_t* pcv, const double* del, const double* azel, double* dant)
{
	double e[3], off[3], cosel = cos(azel[1]);
	int i, j;


	e[0] = sin(azel[0]) * cosel;
	e[1] = cos(azel[0]) * cosel;
	e[2] = sin(azel[1]);

	for (i = 0; i < NFREQ; i++) {
		for (j = 0; j < 3; j++) off[j] = pcv->off[i][j] + del[j];

		dant[i] = -(off[0] * e[0] + off[1] * e[01] + off[2] * e[2]) + interpvar(90.0 - azel[1] * 180.0 / PI, pcv->var[i]);
	}

}

int CBdPPP::sat_yaw(gtime_t time, CString sat, const double* rs, double* exs, double* eys)
{
	double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
	double yaw, cosy, siny, erpv[5] = { 0 };
	int i;

	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

	/* beta and orbit angle */
	for (int i = 0; i < 6; i++)ri[i] = rs[i];

	ri[3] -= OMGe * ri[1];
	ri[4] += OMGe * ri[0];
	cross3(ri, ri + 3, n);
	cross3(rsun, n, p);
	for (int i = 0; i < 3; i++)
	{
		es[i] = rs[i] / sqrt(pow(rs[0], 2) + pow(rs[1], 2) + pow(rs[2], 2));
		esun[i] = rsun[i] / sqrt(pow(rsun[0], 2) + pow(rsun[1], 2) + pow(rsun[2], 2));
		en[i] = n[i] / sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2));
		ep[i] = p[i] / sqrt(pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
	}

	beta = PI / 2.0 - acos(esun[0] * en[0] + esun[1] * en[1] + esun[2] * en[2]);
	E = acos(es[0] * ep[0] + es[1] * ep[1] + es[2] * ep[2]);
	mu = PI / 2.0 + ((es[0] * esun[0] + es[1] * esun[1] + es[2] * esun[2]) <= 0 ? -E : E);
	if (mu < -PI / 2.0) mu += 2.0 * PI;
	else if (mu >= PI / 2.0) mu -= 2.0 * PI;

	/* yaw-angle of satellite */
	if (!yaw_angle( beta, mu, &yaw)) return 0;

	/* satellite fixed x,y-vector */
	cross3(en, es, ex);
	cosy = cos(yaw);
	siny = sin(yaw);
	for (i = 0; i < 3; i++) {
		exs[i] = -siny * en[i] + cosy * ex[i];
		eys[i] = -cosy * en[i] - siny * ex[i];
	}
	return 1;
}

double CBdPPP::tropmapf(gtime_t time, const double pos[], const double azel[], double* mapfw)
{
#ifdef IERS_MODEL
	const double ep[] = { 2000,1,1,12,0,0 };
	double mjd, lat, lon, hgt, zd, gmfh, gmfw;
#endif

	if (pos[2] < -1000.0 || pos[2]>20000.0) {
		if (mapfw) *mapfw = 0.0;
		return 0.0;
	}
#ifdef IERS_MODEL
	mjd = 51544.5 + (timediff(time, epoch2time(ep))) / 86400.0;
	lat = pos[0];
	lon = pos[1];
	hgt = pos[2] - geoidh(pos); /* height in m (mean sea level) */
	zd = PI / 2.0 - azel[1];

	/* call GMF */
	gmf_(&mjd, &lat, &lon, &hgt, &zd, &gmfh, &gmfw);

	if (mapfw) *mapfw = gmfw;
	return gmfh;
#else
	////////////////////////////////////////////////////////PANLIN
	//return nmf(time,pos,azel,mapfw); /* NMF */
	return gmf(time, pos, azel, mapfw); /* GMF */
	////////////////////////////////////////////////////////PANLIN
#endif
}

double CBdPPP::interpvar(double ang, const double* var)
{
	double a = ang / 5.0; /* ang=0-90 */
	int i = (int)a;
	if (i < 0) return var[0]; else if (i >= 18) return var[18];
	return var[i] * (1.0 - a + i) + var[i + 1] * (a - i);
}

double CBdPPP::gmf(gtime_t time, const double pos[], const double azel[], double* gmfw)
{
	const double ep[] = { 2000,1,1,12,0,0 };
	double mjd = 51544.5 + (timediff(gpst2utc(time), epoch2time(ep))) / 86400.0;

	const double ah_mean[55] = {
		+1.2517E+02, +8.503E-01, +6.936E-02, -6.760E+00, +1.771E-01,
		+1.130E-02, +5.963E-01, +1.808E-02, +2.801E-03, -1.414E-03,
		-1.212E+00, +9.300E-02, +3.683E-03, +1.095E-03, +4.671E-05,
		+3.959E-01, -3.867E-02, +5.413E-03, -5.289E-04, +3.229E-04,
		+2.067E-05, +3.000E-01, +2.031E-02, +5.900E-03, +4.573E-04,
		-7.619E-05, +2.327E-06, +3.845E-06, +1.182E-01, +1.158E-02,
		+5.445E-03, +6.219E-05, +4.204E-06, -2.093E-06, +1.540E-07,
		-4.280E-08, -4.751E-01, -3.490E-02, +1.758E-03, +4.019E-04,
		-2.799E-06, -1.287E-06, +5.468E-07, +7.580E-08, -6.300E-09,
		-1.160E-01, +8.301E-03, +8.771E-04, +9.955E-05, -1.718E-06,
		-2.012E-06, +1.170E-08, +1.790E-08, -1.300E-09, +1.000E-10
	};
	const double bh_mean[55] = {
		+0.000E+00, +0.000E+00, +3.249E-02, +0.000E+00, +3.324E-02,
		+1.850E-02, +0.000E+00, -1.115E-01, +2.519E-02, +4.923E-03,
		+0.000E+00, +2.737E-02, +1.595E-02, -7.332E-04, +1.933E-04,
		+0.000E+00, -4.796E-02, +6.381E-03, -1.599E-04, -3.685E-04,
		+1.815E-05, +0.000E+00, +7.033E-02, +2.426E-03, -1.111E-03,
		-1.357E-04, -7.828E-06, +2.547E-06, +0.000E+00, +5.779E-03,
		+3.133E-03, -5.312E-04, -2.028E-05, +2.323E-07, -9.100E-08,
		-1.650E-08, +0.000E+00, +3.688E-02, -8.638E-04, -8.514E-05,
		-2.828E-05, +5.403E-07, +4.390E-07, +1.350E-08, +1.800E-09,
		+0.000E+00, -2.736E-02, -2.977E-04, +8.113E-05, +2.329E-07,
		+8.451E-07, +4.490E-08, -8.100E-09, -1.500E-09, +2.000E-10
	};
	const double ah_amp[55] = {
		-2.738E-01, -2.837E+00, +1.298E-02, -3.588E-01, +2.413E-02,
		+3.427E-02, -7.624E-01, +7.272E-02, +2.160E-02, -3.385E-03,
		+4.424E-01, +3.722E-02, +2.195E-02, -1.503E-03, +2.426E-04,
		+3.013E-01, +5.762E-02, +1.019E-02, -4.476E-04, +6.790E-05,
		+3.227E-05, +3.123E-01, -3.535E-02, +4.840E-03, +3.025E-06,
		-4.363E-05, +2.854E-07, -1.286E-06, -6.725E-01, -3.730E-02,
		+8.964E-04, +1.399E-04, -3.990E-06, +7.431E-06, -2.796E-07,
		-1.601E-07, +4.068E-02, -1.352E-02, +7.282E-04, +9.594E-05,
		+2.070E-06, -9.620E-08, -2.742E-07, -6.370E-08, -6.300E-09,
		+8.625E-02, -5.971E-03, +4.705E-04, +2.335E-05, +4.226E-06,
		+2.475E-07, -8.850E-08, -3.600E-08, -2.900E-09, +0.000E+00
	};
	const double bh_amp[55] = {
		+0.000E+00, +0.000E+00, -1.136E-01, +0.000E+00, -1.868E-01,
		-1.399E-02, +0.000E+00, -1.043E-01, +1.175E-02, -2.240E-03,
		+0.000E+00, -3.222E-02, +1.333E-02, -2.647E-03, -2.316E-05,
		+0.000E+00, +5.339E-02, +1.107E-02, -3.116E-03, -1.079E-04,
		-1.299E-05, +0.000E+00, +4.861E-03, +8.891E-03, -6.448E-04,
		-1.279E-05, +6.358E-06, -1.417E-07, +0.000E+00, +3.041E-02,
		+1.150E-03, -8.743E-04, -2.781E-05, +6.367E-07, -1.140E-08,
		-4.200E-08, +0.000E+00, -2.982E-02, -3.000E-03, +1.394E-05,
		-3.290E-05, -1.705E-07, +7.440E-08, +2.720E-08, -6.600E-09,
		+0.000E+00, +1.236E-02, -9.981E-04, -3.792E-05, -1.355E-05,
		+1.162E-06, -1.789E-07, +1.470E-08, -2.400E-09, -4.000E-10
	};
	const double aw_mean[55] = {
		+5.640E+01, +1.555E+00, -1.011E+00, -3.975E+00, +3.171E-02,
		+1.065E-01, +6.175E-01, +1.376E-01, +4.229E-02, +3.028E-03,
		+1.688E+00, -1.692E-01, +5.478E-02, +2.473E-02, +6.059E-04,
		+2.278E+00, +6.614E-03, -3.505E-04, -6.697E-03, +8.402E-04,
		+7.033E-04, -3.236E+00, +2.184E-01, -4.611E-02, -1.613E-02,
		-1.604E-03, +5.420E-05, +7.922E-05, -2.711E-01, -4.406E-01,
		-3.376E-02, -2.801E-03, -4.090E-04, -2.056E-05, +6.894E-06,
		+2.317E-06, +1.941E+00, -2.562E-01, +1.598E-02, +5.449E-03,
		+3.544E-04, +1.148E-05, +7.503E-06, -5.667E-07, -3.660E-08,
		+8.683E-01, -5.931E-02, -1.864E-03, -1.277E-04, +2.029E-04,
		+1.269E-05, +1.629E-06, +9.660E-08, -1.015E-07, -5.000E-10
	};
	const double bw_mean[55] = {
		+0.000E+00, +0.000E+00, +2.592E-01, +0.000E+00, +2.974E-02,
		-5.471E-01, +0.000E+00, -5.926E-01, -1.030E-01, -1.567E-02,
		+0.000E+00, +1.710E-01, +9.025E-02, +2.689E-02, +2.243E-03,
		+0.000E+00, +3.439E-01, +2.402E-02, +5.410E-03, +1.601E-03,
		+9.669E-05, +0.000E+00, +9.502E-02, -3.063E-02, -1.055E-03,
		-1.067E-04, -1.130E-04, +2.124E-05, +0.000E+00, -3.129E-01,
		+8.463E-03, +2.253E-04, +7.413E-05, -9.376E-05, -1.606E-06,
		+2.060E-06, +0.000E+00, +2.739E-01, +1.167E-03, -2.246E-05,
		-1.287E-04, -2.438E-05, -7.561E-07, +1.158E-06, +4.950E-08,
		+0.000E+00, -1.344E-01, +5.342E-03, +3.775E-04, -6.756E-05,
		-1.686E-06, -1.184E-06, +2.768E-07, +2.730E-08, +5.700E-09
	};
	const double aw_amp[55] = {
		+1.023E-01, -2.695E+00, +3.417E-01, -1.405E-01, +3.175E-01,
		+2.116E-01, +3.536E+00, -1.505E-01, -1.660E-02, +2.967E-02,
		+3.819E-01, -1.695E-01, -7.444E-02, +7.409E-03, -6.262E-03,
		-1.836E+00, -1.759E-02, -6.256E-02, -2.371E-03, +7.947E-04,
		+1.501E-04, -8.603E-01, -1.360E-01, -3.629E-02, -3.706E-03,
		-2.976E-04, +1.857E-05, +3.021E-05, +2.248E+00, -1.178E-01,
		+1.255E-02, +1.134E-03, -2.161E-04, -5.817E-06, +8.836E-07,
		-1.769E-07, +7.313E-01, -1.188E-01, +1.145E-02, +1.011E-03,
		+1.083E-04, +2.570E-06, -2.140E-06, -5.710E-08, +2.000E-08,
		-1.632E+00, -6.948E-03, -3.893E-03, +8.592E-04, +7.577E-05,
		+4.539E-06, -3.852E-07, -2.213E-07, -1.370E-08, +5.800E-09
	};
	const double bw_amp[55] = {
		+0.000E+00, +0.000E+00, -8.865E-02, +0.000E+00, -4.309E-01,
		+6.340E-02, +0.000E+00, +1.162E-01, +6.176E-02, -4.234E-03,
		+0.000E+00, +2.530E-01, +4.017E-02, -6.204E-03, +4.977E-03,
		+0.000E+00, -1.737E-01, -5.638E-03, +1.488E-04, +4.857E-04,
		-1.809E-04, +0.000E+00, -1.514E-01, -1.685E-02, +5.333E-03,
		-7.611E-05, +2.394E-05, +8.195E-06, +0.000E+00, +9.326E-02,
		-1.275E-02, -3.071E-04, +5.374E-05, -3.391E-05, -7.436E-06,
		+6.747E-07, +0.000E+00, -8.637E-02, -3.807E-03, -6.833E-04,
		-3.861E-05, -2.268E-05, +1.454E-06, +3.860E-07, -1.068E-07,
		+0.000E+00, -2.658E-02, -1.947E-03, +7.131E-04, -3.506E-05,
		+1.885E-07, +5.792E-07, +3.990E-08, +2.000E-08, -5.700E-09
	};
	const double aht[] = { 2.53E-5,5.49E-3,1.14E-3 }; /* height correction */

	double dfac[20], P[10][10], aP[55], bP[55], ah[3], aw[3];
	double doy, t, c0h, phh, c11h, c10h, ahm, aha, awm, awa, dm, gmfh;
	double az = azel[0], el = azel[1], lat = pos[0], lon = pos[1], hgt;
	int i, j, k, n, m;
	hgt = pos[2];

	doy = mjd - 44239.0 + 1 - 28;
	/* parameter t */
	t = sin(lat);
	n = 9; m = 9; /* degree n and order m EGM */

	/* determine n! (faktorielle)  moved by 1 */
	dfac[0] = 1;
	for (i = 1; i <= 2 * n + 1; i++) dfac[i] = dfac[i - 1] * (i);

	Pnm(9, 9, sin(lat), P[0]);

	/* spherical harmonics */
	for (i = 0, k = 0; i <= 9; i++) {
		for (j = 0; j <= i; j++, k++) {
			aP[k] = P[i][j] * cos(j * lon);
			bP[k] = P[i][j] * sin(j * lon);
		}
	}

	/****************** hydrostatic mapping function *****************/
	ah[1] = 0.0029; c0h = 0.062;
	if (lat < 0) { /* southern hemisphere */
		phh = PI; c11h = 0.007; c10h = 0.002;
	}
	else {       /* northern hemisphere */
		phh = 0;  c11h = 0.005; c10h = 0.001;
	}
	ah[2] = c0h + ((cos(doy / 365.25 * 2 * PI + phh) + 1) * c11h / 2 + c10h) * (1 - cos(lat));
	ahm = 0.0; aha = 0.0;
	for (i = 0; i < 55; i++) {
		ahm += (ah_mean[i] * aP[i] + bh_mean[i] * bP[i]) * 1E-5;
		aha += (ah_amp[i] * aP[i] + bh_amp[i] * bP[i]) * 1E-5;
	}
	ah[0] = ahm + aha * cos(doy / 365.25 * 2.0 * PI);
	dm = (1.0 / sin(el) - mapf(el, aht[0], aht[1], aht[2])) * (hgt - geoidh(pos)) / 1E3;
	gmfh = mapf(el, ah[0], ah[1], ah[2]) + dm;

	/********************* wet mapping function *********************/
	aw[1] = 0.00146; aw[2] = 0.04391;
	awm = 0.0; awa = 0.0;
	for (i = 0; i < 55; i++) {
		awm += (aw_mean[i] * aP[i] + bw_mean[i] * bP[i]) * 1E-5;
		awa += (aw_amp[i] * aP[i] + bw_amp[i] * bP[i]) * 1E-5;
	}
	aw[0] = awm + awa * cos(doy / 365.25 * 2.0 * PI);
	if (gmfw) *gmfw = mapf(el, aw[0], aw[1], aw[2]);

	return gmfh;
}

int CBdPPP::yaw_angle(double beta, double mu, double* yaw)
{
	*yaw = yaw_nominal(beta, mu);
	return 1;
}

double CBdPPP::geoidh(const double* pos)
{
	double posd[2], h;

	posd[1] = pos[1] * 180.0/PI; posd[0] = pos[0] * 180.0 / PI; if (posd[1] < 0.0) posd[1] += 360.0;

	h = geoidh_emb(posd);

	return h;
}

double CBdPPP::geoidh_emb(const double* pos)
{
	const double dlon = 1.0, dlat = 1.0;
	double a, b, y[4];
	int i1, i2, j1, j2;

	a = (pos[1] - range[0]) / dlon;
	b = (pos[0] - range[2]) / dlat;
	i1 = (int)a; a -= i1; i2 = i1 < 360 ? i1 + 1 : i1;
	j1 = (int)b; b -= j1; j2 = j1 < 180 ? j1 + 1 : j1;
	y[0] = geoid[i1][j1];
	y[1] = geoid[i2][j1];
	y[2] = geoid[i1][j2];
	y[3] = geoid[i2][j2];
	return interpb(y, a, b);
}

void CBdPPP::Pnm(int n, int m, double x, double* P)
{
	int i = 0, j = 0, k = 0;
	double ir = 0, sum = 0, * dfac;

	/* determine n! (faktorielle)  moved by 1 */
	dfac = (double*)malloc(sizeof(double) * (2 * n + 2));
	dfac[0] = 1;
	for (i = 1; i <= 2 * n + 1; i++) dfac[i] = dfac[i - 1] * (i);

	/* determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62) */
	for (i = 0; i <= n; i++) {
		for (j = 0; j <= min(i, m); j++) {
			ir = (int)((i - j) / 2);
			sum = 0;
			for (k = 0; k <= ir; k++) {
				sum += ((int)pow((double)(-1.0), k)) * dfac[2 * i - 2 * k] / dfac[k] / dfac[i - k] / dfac[i - j - 2 * k] * pow(x, i - j - 2 * k);
			}
			/*  Legendre functions moved by 1 */
			P[i * (m + 1) + j] = 1.0 / pow((double)(2.0), i) * sqrt(pow(1 - x * x, j)) * sum;
		}
	}

	free(dfac);
}

double CBdPPP::mapf(double el, double a, double b, double c)
{
	double sinel = sin(el);
	return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))));
}

double CBdPPP::interpb(const double* y, double a, double b)
{
	return y[0] * (1.0 - a) * (1.0 - b) + y[1] * a * (1.0 - b) + y[2] * (1.0 - a) * b + y[3] * a * b;
}

 void CBdPPP::cross3(const double* a, const double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

double CBdPPP::yaw_nominal(double beta, double mu)
 {
	 if (fabs(beta) < 1E-12 && fabs(mu) < 1E-12) return PI;
	 return atan2(-tan(beta), sin(mu)) + PI;
 }

double CBdPPP::varerr(CString sat, double el, int freq, int type, const double* pos)   ////////////////////////////PANLIN  add parameter "const double *pos" for uncombined PPP with respect to BDS 
{
	double fact = 1.0, sinel = sin(el);
	int prn=_ttoi(sat.Right(2)); /////////PANLIN

	///////////
	if (type == 1) fact *= eratio[freq == 0 ? 0 : 1];


	//////////////Post Processing

	////////////IFLC model
	fact *= 3.0;

	///////////////////////////////////////////////////////////////////PANLIN
	///*
	fact *= 2.0;
	if ((prn <= 5) || (prn == 59) || (prn == 60) || (prn == 61))   //GEO Sat
	{
		fact *= 2.0;
	}
	//if ( (opt->mode==PMODE_PPP_KINEMA) && (sys==SYS_CMP) )  fact*=2.0;
	//if (opt->ionoopt==IONOOPT_EST)
	//{
	//	/*
	//	if (sys==SYS_GAL) fact*=20.0;
	//	if (sys==SYS_CMP)
	//	{
	//		if ( (pos[0]*R2D>-55) && (pos[0]*R2D<55) && (pos[1]*R2D>70) && (pos[1]*R2D<150) ) fact*=1.0;
	//		else fact*=20.0;
	//	}
	//	*/
	//	if ( (sys==SYS_CMP) && (type==1) ) fact*=100.0;
	//}
	//*/
	///////////////////////////////////////////////////////////////////PANLIN

	return SQR(fact * err[1]) + SQR(fact *err[2] / sinel);
}


double CBdPPP::sqrtvar(double var)
{
	return var > 0 ? sqrt(var) : -sqrt(-var);
}




