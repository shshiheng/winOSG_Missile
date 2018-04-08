#include "Earth.h"
#include "EarthConstant.h"
#include "MathConstant.h"

#include <math.h>
#include "UserMechanics.h"
#include "TimeConstant.h"

void CEarth::Setup(const double MechanicsTime)
{
	this->JulianCentury = MechanicsTime;
	this->SetRei();
}
void CEarth::Step(const double TimeStep)
{
	this->JulianCentury += TimeStep/86400.0/36525.0;
	this->SetRei();
}

void CEarth::SetRei()
{
	//岁差矩阵
	this->PrecessionZeta = 2306.2181/3600/180*PI*this->JulianCentury + 0.30188/3600/180*PI*this->JulianCentury*this->JulianCentury + 0.017998/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	this->PrecessionZi = 2306.2181/3600/180*PI*this->JulianCentury + 1.09468/3600/180*PI*this->JulianCentury*this->JulianCentury + 0.018203/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	this->PrecessionThet = 2004.3109/3600/180*PI*this->JulianCentury - 0.42665/3600/180*PI*this->JulianCentury*this->JulianCentury - 0.041833/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;

	double Rz_zi[3][3],Ry_thet[3][3],Rz_zeta[3][3],temp[3][3],PM[3][3];
	GetRz(-this->PrecessionZi,Rz_zi);
	GetRy(this->PrecessionThet,Ry_thet);
	GetRz(-this->PrecessionZeta,Rz_zeta);

	MatrixProduct(3, 3, 3, &Rz_zi[0][0], &Ry_thet[0][0], &temp[0][0]);
	MatrixProduct(3, 3, 3, &temp[0][0], &Rz_zeta[0][0], &PM[0][0]);

	//章动矩阵
	double Omiga,F,D,M,det_eps;
	Omiga = (125+2.0/60+40.280/3600)/180*PI - (1934+8.0/60+10.539/3600)/180*PI*this->JulianCentury + 7.455/3600/180*PI*this->JulianCentury*this->JulianCentury + 0.008/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	F = (93+16.0/60+18.877/3600)/180*PI + (483202+1.0/60+3.137/3600)/180*PI*this->JulianCentury - 13.257/3600/180*PI*this->JulianCentury*this->JulianCentury + 0.011/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	D = (297+51.0/60+1.307/3600)/180*PI + (445267+6.0/60+41.328/3600)/180*PI*this->JulianCentury - 6.891/3600/180*PI*this->JulianCentury*this->JulianCentury + 0.019/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	M = (357+31.0/60+39.804/3600)/180*PI + (35999+3.0/60+1.224/3600)/180*PI*this->JulianCentury - 0.577/3600/180*PI*this->JulianCentury*this->JulianCentury - 0.012/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;

	this->NutationInLongitude = -(17.1996/3600/180*PI + 0.01742/3600/180*PI*this->JulianCentury)*sin(Omiga)
		+ (0.2062/3600/180*PI + 0.00002/3600/180*PI*this->JulianCentury)*sin(2*Omiga)
		- (1.3187/3600/180*PI + 0.00016/3600/180*PI*this->JulianCentury)*sin(2*F-2*D+2*Omiga)
		+ (0.1426/3600/180*PI - 0.00034/3600/180*PI*this->JulianCentury)*sin(M)
		- (0.2274/3600/180 + 0.00002*this->JulianCentury)*sin(2*F+2*Omiga);
	det_eps = (9.2025/3600/180*PI + 0.00089/3600/180*PI*this->JulianCentury)*cos(Omiga)
		- (0.0895/3600/180*PI - 0.00005/3600/180*PI*this->JulianCentury)*cos(2*Omiga)
		+ (0.5736/3600/180*PI - 0.00031/3600/180*PI*this->JulianCentury)*cos(2*F-2*D+2*Omiga)
		+ (0.0977/3600/180*PI - 0.00005/3600/180*PI*this->JulianCentury)*cos(2*F+2*Omiga)
		+ (0.0054/3600/180*PI - 0.00001/3600/180*PI*this->JulianCentury)*cos(M);
	this->MeanObliquityEcliptic = (23+26.0/60+21.448/3600)/180*PI - 46.8150/3600/180*PI*this->JulianCentury - 0.00059/3600/180*PI*this->JulianCentury*this->JulianCentury + 0.001813/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	this->TrueObliquityEcliptic = this->MeanObliquityEcliptic + det_eps;

	double Rx_eps[3][3],Rz_det_pusai[3][3],Rx_eps_A[3][3],NM[3][3];
	GetRx(-this->TrueObliquityEcliptic,Rx_eps);
	GetRz(-this->NutationInLongitude,Rz_det_pusai);
	GetRx(this->MeanObliquityEcliptic,Rx_eps_A);
	MatrixProduct(3, 3, 3, &Rx_eps[0][0], &Rz_det_pusai[0][0], &temp[0][0]);
	MatrixProduct(3, 3, 3, &temp[0][0], &Rx_eps_A[0][0], &NM[0][0]);

	//格林威治角
	double EAR,EE;
	EAR = 2*PI*(0.7790572732640l + 1.00273781191135448l*(this->JulianCentury*36525- UTC_to_TDT/86400));
	EE = this->NutationInLongitude*cos(this->MeanObliquityEcliptic) + 0.00264/3600/180*PI*sin(Omiga) + 0.000063/3600/180*PI*sin(2*Omiga);
	this->GreenwichAngle = EAR + EE + 0.014506/3600/180*PI + 4612.15739966l/3600/180*PI*this->JulianCentury + 1.39667721l/3600/180*PI*this->JulianCentury*this->JulianCentury - 0.00009344l/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury + 0.00001882l/3600/180*PI*this->JulianCentury*this->JulianCentury*this->JulianCentury*this->JulianCentury;
	while(this->GreenwichAngle > 2*PI)		this->GreenwichAngle = this->GreenwichAngle - 2*PI;
	double SM[3][3];
	GetRz(this->GreenwichAngle,SM);

	//极移，2012年1月的极移数据，xp = 100,yp = 250,单位milliarcsec
	double xp = 100.0/1000/3600/180*PI;
	double yp = 250.0/1000/3600/180*PI;
	double Ry_xp[3][3],Rx_yp[3][3],WM[3][3];
	GetRy(-xp,Ry_xp);
	GetRx(-yp,Rx_yp);
	MatrixProduct(3, 3, 3, &Ry_xp[0][0], &Rx_yp[0][0], &WM[0][0]);

	//Rei = WM*SM*NM*PM
	double temp1[3][3];
	MatrixProduct(3, 3, 3, &WM[0][0], &SM[0][0], &temp[0][0]);
	MatrixProduct(3, 3, 3, &temp[0][0], &NM[0][0], &temp1[0][0]);
	MatrixProduct(3, 3, 3, &temp1[0][0], &PM[0][0], &this->Rei[0][0]);
}