
// BeiDouPPPDlg.cpp: 实现文件
//

#include "pch.h"
#include "framework.h"
#include "BeiDouPPP.h"
#include "BeiDouPPPDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CBeiDouPPPDlg 对话框



CBeiDouPPPDlg::CBeiDouPPPDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_BEIDOUPPP_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CBeiDouPPPDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CBeiDouPPPDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_ReadSP3, &CBeiDouPPPDlg::OnBnClickedReadsp3)
	ON_BN_CLICKED(IDC_ReadEph, &CBeiDouPPPDlg::OnBnClickedReadeph)
	ON_BN_CLICKED(IDC_ReadObs, &CBeiDouPPPDlg::OnBnClickedReadobs)
	ON_BN_CLICKED(IDC_ReadClk, &CBeiDouPPPDlg::OnBnClickedReadclk)
	ON_BN_CLICKED(IDC_ReadAtx, &CBeiDouPPPDlg::OnBnClickedReadatx)
	ON_BN_CLICKED(IDC_ReadOtl, &CBeiDouPPPDlg::OnBnClickedReadotl)
	ON_BN_CLICKED(IDC_ReadErp, &CBeiDouPPPDlg::OnBnClickedReaderp)
	ON_BN_CLICKED(IDC_ReadBsx, &CBeiDouPPPDlg::OnBnClickedReadbsx)
	ON_BN_CLICKED(IDC_MainProcess, &CBeiDouPPPDlg::OnBnClickedMainprocess)
	ON_BN_CLICKED(IDC_OutSol, &CBeiDouPPPDlg::OnBnClickedOutsol)
END_MESSAGE_MAP()


// CBeiDouPPPDlg 消息处理程序

BOOL CBeiDouPPPDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CBeiDouPPPDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CBeiDouPPPDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CBeiDouPPPDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CBeiDouPPPDlg::OnBnClickedReadeph()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22p"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22p"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/BRDM00DLR_S_20220320000_01D_MN.22p");
	m_CBdPPP.ReadNavData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReadobs()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/jfng0320.22o");
	m_CBdPPP.ReadObsData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReadsp3()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/WUM0MGXFIN_20220320000_01D_15M_ORB.SP3");
	m_CBdPPP.ReadPreEphData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReadclk()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/WUM0MGXFIN_20220320000_01D_30S_CLK.CLK");
	m_CBdPPP.ReadPreClkData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReadatx()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/igs20_2274.atx");
	m_CBdPPP.ReadPcvData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReadotl()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/FES2004.BLQ");
	m_CBdPPP.ReadOtlData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReaderp()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/WUM0MGXFIN_20220320000_01D_01D_ERP.ERP");
	m_CBdPPP.ReadErpData(strFileName);
}

void CBeiDouPPPDlg::OnBnClickedReadbsx()
{
	/*
	CFileDialog dlgFile(TRUE, _T("22o"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.22o"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	//NavFile = strFileName;
	CString strFileName("..//..//PPP_Data/CAS0MGXRAP_20220320000_01D_01D_DCB.BSX");
	m_CBdPPP.ReadDcbData(strFileName);
}


void CBeiDouPPPDlg::OnBnClickedMainprocess()
{
	/*
	CFileDialog dlgFile(TRUE, _T("dat"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.dat"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	m_CBdPPP.ReadNavData(_T("..//..//PPP_Data/BRDM00DLR_S_20220320000_01D_MN.22p"));

	m_CBdPPP.ReadObsData(_T("..//..//PPP_Data/jfng0320.22o"));

	m_CBdPPP.ReadPreEphData(_T("..//..//PPP_Data/WUM0MGXFIN_20220320000_01D_15M_ORB.SP3"));

	m_CBdPPP.ReadPreClkData(_T("..//..//PPP_Data/WUM0MGXFIN_20220320000_01D_30S_CLK.CLK"));

	m_CBdPPP.ReadPcvData(_T("..//..//PPP_Data/igs20_2274.atx"));

	m_CBdPPP.ReadOtlData(_T("..//..//PPP_Data/FES2004.BLQ"));

	m_CBdPPP.ReadErpData(_T("..//..//PPP_Data/WUM0MGXFIN_20220320000_01D_01D_ERP.ERP"));

	m_CBdPPP.ReadDcbData(_T("..//..//PPP_Data/CAS0MGXRAP_20220320000_01D_01D_DCB.BSX"));

	m_CBdPPP.pppos(_T("..//..//PPP_Data/testMat.dat"));
}


void CBeiDouPPPDlg::OnBnClickedOutsol()
{
	/*
	CFileDialog dlgFile(TRUE, _T("pos"), NULL,
		OFN_ALLOWMULTISELECT | OFN_EXPLORER,
		_T("(文本文件)|*.pos"));
	if (dlgFile.DoModal() == IDCANCEL) return;
	CString strFileName = dlgFile.GetPathName();
	*/
	CString strFileName("..//..//PPP_Data/myPos.pos");
	m_CBdPPP.OutputSol(strFileName);
}
