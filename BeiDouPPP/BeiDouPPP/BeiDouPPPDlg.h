
// BeiDouPPPDlg.h: 头文件
//

#pragma once
#include "BdPPP.h"

// CBeiDouPPPDlg 对话框
class CBeiDouPPPDlg : public CDialogEx
{
// 构造
public:
	CBeiDouPPPDlg(CWnd* pParent = nullptr);	// 标准构造函数

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_BEIDOUPPP_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;
	CBdPPP m_CBdPPP;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedReadsp3();
	afx_msg void OnBnClickedReadeph();
	afx_msg void OnBnClickedReadobs();
	afx_msg void OnBnClickedReadclk();
	afx_msg void OnBnClickedReadatx();
	afx_msg void OnBnClickedReadotl();
	afx_msg void OnBnClickedReaderp();
	afx_msg void OnBnClickedReadbsx();
	afx_msg void OnBnClickedMainprocess();
	afx_msg void OnBnClickedOutsol();
};
