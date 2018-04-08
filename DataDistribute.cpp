//File Satellite contains the realization of Satellite functions
//File created @ 2009-11-19
#include <stdio.h>
#include <Math.h>

#include "DataDistribute.h"
#include "Hud.h"
#include "HLBToFixXYZ.h"
#include "Satellite.h"
#include "UserMath.h"
#include "CreateDevice.h"
#include "DebugDefine.h"
#include "RTI_NET.h"
#include "RTCHeadDefine.h"
#include "FederalNo.h"
#include "UnitConstant.h"
#include "disEarth.h"
#include "Earth.h"
#include "disMotion.h"
#include "UnitConstant.h"

#ifdef CONNECT_ARCHIVE
#include "ArchiveBUSChannelNo.h"
#endif

#ifdef USE_AS_OSG
#include "OSGBUSChannelNo.h"
#endif

#define MAX_MAINBODY_NUMBER 5
#define MAX_DATATYPE_NUMBER 50
static double Arithmometer_RTI=0;
MainBody *pCurrentSatellite = NULL;
extern MainBody *pMainBodyInstance[MAX_MAINBODY_NUMBER];
SceneManipulator *pCurrentManipulator = NULL;
Hand *pCurrentHand = NULL;
static CMotionAgent		clsTimer;
static CEarth			clsEarth;

double SimulationTime_9A = 0;
double SimulationTime_9B = 0;
double Satellite9APosition_X = 0;
double Satellite9APosition_Y = 0;
double Satellite9APosition_Z = 0;
double Satellite9AQuaternionAttitude_Q0 = 0;
double Satellite9AQuaternionAttitude_Q1 = 0;
double Satellite9AQuaternionAttitude_Q2 = 0;
double Satellite9AQuaternionAttitude_Q3 = 0;
double Satellite9ASpeed[3] = {0.0};
double Satellite9AAngleSpeed[3] = {0.0};
double Satellite9AEular[3] = {0.0};
double Satellite9BPosition_X = 0;
double Satellite9BPosition_Y = 0;
double Satellite9BPosition_Z = 0;
double Satellite9BQuaternionAttitude_Q0 = 0;
double Satellite9BQuaternionAttitude_Q1 = 0;
double Satellite9BQuaternionAttitude_Q2 = 0;
double Satellite9BQuaternionAttitude_Q3 = 0;
double Satellite9BSpeed[3] = {0.0};
double Satellite9BAngleSpeed[3] = {0.0};
double Satellite9BEular[3] = {0.0};	
double SatelliteRelativePosition = 0.0;




void RegisterManipulatorUpdate(SceneManipulator *pManipulator)
{
	pCurrentManipulator = pManipulator ;
}

void RegisterHandUpdate(Hand *pHand)
{
	pCurrentHand = pHand;
}

void RegisterSatelliteUpdate(MainBody *pSatellite)
{
	pCurrentSatellite = pSatellite;
}


bool Satellite_A_BodyShow = true;
bool Satellite_B_BodyShow = true;

void CollectSatelliteData(void)	
{
	static double SatellitePosition_X = 0;
	static double SatellitePosition_Y = 0;
	static double SatellitePosition_Z = 0;
	static double SatellitePosition_Fix[3] = {0.0};
	static double SatellitePosition_Initial[3] = {0.0};
	static double SatelliteSpeed[3] = {0.0};
	static double SatelliteAngleSpeed[3] = {0.0};
	static double SatelliteEular[3] = {0.0};
	static double SatelliteQuaternionAttitude_Q0 = 0;
	static double SatelliteQuaternionAttitude_Q1 = 0;
	static double SatelliteQuaternionAttitude_Q2 = 0;
	static double SatelliteQuaternionAttitude_Q3 = 0;
	static double mTranMatrix[3][3] = {0.0};
	static double Q[4] = {0.0};
	static double Relative_R1[3] = {0.0};
	static double Relative_R[3] = {0.0};
	static double Satellite9A_R[3] = {0.0};
	static double Satellite9A_V[3] = {0.0};
	static long lEarthSensorViewer = 0;
	static long lEarthSensorDeviceViewer[2] = {0};
	static long lSunSensorViewer = 0;
	static long lSunSensorDeviceViewer[2] = {0};
	double GyroValue = 0.0;
	static double Satellite_A_GyroDeviceValue[6] = {0.0};
	static double Satellite_B_GyroDeviceValue[6] = {0.0};
	double EarthSensorValues[2] = {0.0};
	static double EarthSensorRollValues[2] = {0.0};
	static double EarthSensorPitchValues[2] = {0.0};
	static double SunSensorValues[2] = {0.0};
	static double Satellite_A_SunSensorValues[2] = {0.0};
	static double Satellite_B_SunSensorValues[2] = {0.0};
	static double MagnetometerValue[3] = {0.0};
	static double Satellite_A_MagnetometerValue[3] = {0.0};
	static double Satellite_B_MagnetometerValue[3] = {0.0};
	static double ThrusterValue = 0.0;
	static double WheelValues[2] = {0.0};
	static double Satellite_A_WheelSpeed[4] = {0.0};
	static double Satellite_A_WheelAccelerate[4] = {0.0};
	static double Satellite_B_WheelSpeed[4] = {0.0};
	static double Satellite_B_WheelAccelerate[4] = {0.0};
	static double CoilValue[3] = {0.0};
	static double Satellite_A_CoilValue[3] = {0.0};
	static double Satellite_B_CoilValue[3] = {0.0};
	static double q[4] = {0.0};
	static double Satellite_A_Q0[4] = {0.0};
	static double Satellite_A_Q1[4] = {0.0};
	static double Satellite_B_Q0[4] = {0.0};
	static double Satellite_B_Q1[4] = {0.0};
	static SGeographyElements	HLB;
	/*RTI*/
	unsigned char j = 0;
	unsigned char ucTotalNo = 0;
	struct SRTIMessage stcRTIBuffers[RTI_BUFFER_SIZE];

	CEarthAgent clsEarthSatellite; 

	ucTotalNo = CollectRTIMessages(ulDBArchiveChannelNoforArchiveAgent,stcRTIBuffers);
	if ( ucTotalNo < 1 ) 
	{
		return;
	}
	for(j = 0;j<ucTotalNo;j++)
	{
		char chrDataType[MAX_DATATYPE_NUMBER] = "\0";
		long lCraftNo = 0;
		long lDeviceNo = 0;
		double SimulationTime = 0;
		unsigned char ucLen = 0;
		unsigned char ucCount = 0;
		unsigned char ucLength = 0;
		static unsigned char uc9ABoxNoInOSG = 0;
		static unsigned char uc9BBoxNoInOSG = 0;
		static unsigned char ucFirstTime = 1;
		std::string strCompare;

		ucLength += 5;
		memcpy(&ucLen,&stcRTIBuffers[j].chrMessage[ucLength],1);										ucLength += 1;
		memcpy(&chrDataType,&stcRTIBuffers[j].chrMessage[ucLength],ucLen);								ucLength += ucLen;
		if(strcmp(chrDataType,"MissileTrajectory") != 0 && strcmp(chrDataType,"MissileAttitudes") != 0 && strcmp(chrDataType,"MissileTrajectoryToOSG") != 0 && strcmp(chrDataType,"MissileAttitudesToOSG") != 0
			 && strcmp(chrDataType,"DecoyTrajectory") != 0 && strcmp(chrDataType,"DecoyAttitudes") != 0 && strcmp(chrDataType,"DecoyTrajectoryToOSG") != 0 && strcmp(chrDataType,"DecoyAttitudesToOSG") != 0)
			continue;
		if(strcmp(chrDataType,"MissileTrajectory") == 0 || strcmp(chrDataType,"MissileTrajectoryToOSG") == 0  || strcmp(chrDataType,"DecoyTrajectory") == 0 || strcmp(chrDataType,"DecoyTrajectoryToOSG") == 0)
		{
			memcpy(&lCraftNo,&stcRTIBuffers[j].chrMessage[ucLength],4);									ucLength += 4;			
			memcpy(&SimulationTime,&stcRTIBuffers[j].chrMessage[ucLength],8);							ucLength += 8;
			memcpy((void*)&SatellitePosition_X,&stcRTIBuffers[j].chrMessage[ucLength],8);				ucLength += 8;
			memcpy((void*)&SatellitePosition_Y,&stcRTIBuffers[j].chrMessage[ucLength],8);				ucLength += 8;
			memcpy((void*)&SatellitePosition_Z,&stcRTIBuffers[j].chrMessage[ucLength],8);				ucLength += 8;
			if (strcmp(chrDataType,"DecoyTrajectory") == 0)												ucLength += 24;

			memcpy((void*)&SatelliteSpeed[0],&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;
			memcpy((void*)&SatelliteSpeed[1],&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;
			memcpy((void*)&SatelliteSpeed[2],&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;
		}
		else if(strcmp(chrDataType,"MissileAttitudes") == 0 || strcmp(chrDataType,"MissileAttitudesToOSG") == 0 || strcmp(chrDataType,"DecoyAttitudes") == 0 || strcmp(chrDataType,"DecoyAttitudesToOSG") == 0)
		{
			memcpy((void*)&lCraftNo,&stcRTIBuffers[j].chrMessage[ucLength],4);							ucLength += 4;			
			memcpy((void*)&SimulationTime,&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;
			memcpy((void*)&SatelliteQuaternionAttitude_Q0,&stcRTIBuffers[j].chrMessage[ucLength],8);	ucLength += 8;
			memcpy((void*)&SatelliteQuaternionAttitude_Q1,&stcRTIBuffers[j].chrMessage[ucLength],8);	ucLength += 8;
			memcpy((void*)&SatelliteQuaternionAttitude_Q2,&stcRTIBuffers[j].chrMessage[ucLength],8);	ucLength += 8;
			memcpy((void*)&SatelliteQuaternionAttitude_Q3,&stcRTIBuffers[j].chrMessage[ucLength],8);	ucLength += 8;

			memcpy((void*)&SatelliteEular[0],&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;
			memcpy((void*)&SatelliteEular[1],&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;
			memcpy((void*)&SatelliteEular[2],&stcRTIBuffers[j].chrMessage[ucLength],8);					ucLength += 8;

			memcpy((void*)&SatelliteAngleSpeed[0],&stcRTIBuffers[j].chrMessage[ucLength],8);			ucLength += 8;
			memcpy((void*)&SatelliteAngleSpeed[1],&stcRTIBuffers[j].chrMessage[ucLength],8);			ucLength += 8;
			memcpy((void*)&SatelliteAngleSpeed[2],&stcRTIBuffers[j].chrMessage[ucLength],8);			ucLength += 8;
		}

	/*导弹轨道姿态参数*/
		if((strcmp(chrDataType,"MissileAttitudesToOSG")==0)&&lCraftNo == 0)
		{
			SimulationTime_9A = SimulationTime;
			Satellite9AQuaternionAttitude_Q0 = SatelliteQuaternionAttitude_Q1;
			Satellite9AQuaternionAttitude_Q1 = SatelliteQuaternionAttitude_Q2;
			Satellite9AQuaternionAttitude_Q2 = SatelliteQuaternionAttitude_Q3;
			Satellite9AQuaternionAttitude_Q3 = SatelliteQuaternionAttitude_Q0;
			Satellite9AEular[0] = SatelliteEular[0];
			Satellite9AEular[1] = SatelliteEular[1];
			Satellite9AEular[2] = SatelliteEular[2];
			Satellite9AAngleSpeed[0] = SatelliteAngleSpeed[0];
			Satellite9AAngleSpeed[1] = SatelliteAngleSpeed[1];
			Satellite9AAngleSpeed[2] = SatelliteAngleSpeed[2];
		}	
		else if((strcmp(chrDataType,"MissileTrajectoryToOSG")==0)&&lCraftNo == 0)
		{			
			Satellite9APosition_X = SatellitePosition_X;
			Satellite9APosition_Y = SatellitePosition_Y;		
			Satellite9APosition_Z = SatellitePosition_Z;

			Satellite9ASpeed[0] = SatelliteSpeed[0];
			Satellite9ASpeed[1] = SatelliteSpeed[1];
			Satellite9ASpeed[2] = SatelliteSpeed[2];
		}
		
		if(strcmp(chrDataType,"MissileAttitudesToOSG")==0)
		{
			pCurrentSatellite->clsActiveDevices[0].vPosition._v[0] = Satellite9APosition_X;
			pCurrentSatellite->clsActiveDevices[0].vPosition._v[1] = Satellite9APosition_Y;
			pCurrentSatellite->clsActiveDevices[0].vPosition._v[2] = Satellite9APosition_Z;
			pCurrentSatellite->clsActiveDevices[0].vSpeed = osg::Vec3(Satellite9ASpeed[0],Satellite9ASpeed[1],Satellite9ASpeed[2]);

			pCurrentSatellite->clsActiveDevices[0].qQuaternionAttitude._v[0] = Satellite9AQuaternionAttitude_Q0;
			pCurrentSatellite->clsActiveDevices[0].qQuaternionAttitude._v[1] = Satellite9AQuaternionAttitude_Q1;
			pCurrentSatellite->clsActiveDevices[0].qQuaternionAttitude._v[2] = Satellite9AQuaternionAttitude_Q2;
			pCurrentSatellite->clsActiveDevices[0].qQuaternionAttitude._v[3] = Satellite9AQuaternionAttitude_Q3;

			pCurrentSatellite->clsActiveDevices[0].vEulerAttitude._v[0] = Satellite9AEular[0];
			pCurrentSatellite->clsActiveDevices[0].vEulerAttitude._v[1] = Satellite9AEular[1];
			pCurrentSatellite->clsActiveDevices[0].vEulerAttitude._v[2] = Satellite9AEular[2];

			pCurrentSatellite->clsActiveDevices[0].vAngularSpeed._v[0] = Satellite9AAngleSpeed[0];
			pCurrentSatellite->clsActiveDevices[0].vAngularSpeed._v[1] = Satellite9AAngleSpeed[1];
			pCurrentSatellite->clsActiveDevices[0].vAngularSpeed._v[2] = Satellite9AAngleSpeed[2];

			OrbitUpdate(pCurrentSatellite->clsActiveDevices[0].vPosition,pCurrentSatellite->clsActiveDevices[80].pNodeModel->asGeode()->getDrawable(0));
		}		

		if(strcmp(chrDataType,"MissileAttitudes")==0)
		{
			{
				char chrTemp[100]={'\0'};
				sprintf(chrTemp,"Missile Position:\nRx(m):%lf\nRy(m):%lf\nRz(m):%lf\n",SatellitePosition_X,SatellitePosition_Y,SatellitePosition_Z);
				WriteToHud(0,0,chrTemp);
			}
			{
				char chrTemp[100]={'\0'};
				sprintf(chrTemp,"Missile Speed:\nVx(m/s):%lf\nVy(m/s):%lf\nVz(m/s):%lf\n", SatelliteSpeed[0], SatelliteSpeed[1], SatelliteSpeed[2]);
				WriteToHud(0,1,chrTemp);
			}
			{
				char chrTemp[100]={'\0'};
				sprintf(chrTemp,"EularAngle:\nRoll(deg) :%lf\nPitch(deg):%lf\nYew(deg)  :%lf\n",SatelliteEular[0],SatelliteEular[1],SatelliteEular[2]);
				WriteToHud(0,2,chrTemp);
			}
			{	
				char chrTemp[100]={'\0'};
				sprintf(chrTemp,"AngleSpeed:\nX(deg/s):%lf\nY(deg/s):%lf\nZ(deg/s):%lf\n",SatelliteAngleSpeed[0],SatelliteAngleSpeed[1],SatelliteAngleSpeed[2]);
				WriteToHud(0,3,chrTemp);
			}
		}

		for (unsigned i = 0; i < 10; i++)
		{
			if((strcmp(chrDataType,"DecoyAttitudesToOSG")==0) && lCraftNo == i)
			{
				SimulationTime_9B = SimulationTime;
				Satellite9BQuaternionAttitude_Q0 = SatelliteQuaternionAttitude_Q1;
				Satellite9BQuaternionAttitude_Q1 = SatelliteQuaternionAttitude_Q2;
				Satellite9BQuaternionAttitude_Q2 = SatelliteQuaternionAttitude_Q3;
				Satellite9BQuaternionAttitude_Q3 = SatelliteQuaternionAttitude_Q0;
				Satellite9BEular[0] = SatelliteEular[0];
				Satellite9BEular[1] = SatelliteEular[1];
				Satellite9BEular[2] = SatelliteEular[2];
				Satellite9BAngleSpeed[0] = SatelliteAngleSpeed[0];
				Satellite9BAngleSpeed[1] = SatelliteAngleSpeed[1];
				Satellite9BAngleSpeed[2] = SatelliteAngleSpeed[2];
			}	
			else if((strcmp(chrDataType,"DecoyTrajectoryToOSG")==0) && lCraftNo == i)
			{			
				Satellite9BPosition_X = SatellitePosition_X;
				Satellite9BPosition_Y = SatellitePosition_Y;		
				Satellite9BPosition_Z = SatellitePosition_Z;

				Satellite9BSpeed[0] = SatelliteSpeed[0];
				Satellite9BSpeed[1] = SatelliteSpeed[1];
				Satellite9BSpeed[2] = SatelliteSpeed[2];
			}

			if(strcmp(chrDataType,"DecoyAttitudesToOSG")==0 && lCraftNo == i)
			{
				pCurrentSatellite->clsActiveDevices[i+1].vPosition._v[0] = Satellite9BPosition_X;
				pCurrentSatellite->clsActiveDevices[i+1].vPosition._v[1] = Satellite9BPosition_Y;
				pCurrentSatellite->clsActiveDevices[i+1].vPosition._v[2] = Satellite9BPosition_Z;
				pCurrentSatellite->clsActiveDevices[i+1].vSpeed = osg::Vec3(Satellite9BSpeed[0],Satellite9BSpeed[1],Satellite9BSpeed[2]);

				pCurrentSatellite->clsActiveDevices[i+1].qQuaternionAttitude._v[0] = Satellite9BQuaternionAttitude_Q0;
				pCurrentSatellite->clsActiveDevices[i+1].qQuaternionAttitude._v[1] = Satellite9BQuaternionAttitude_Q1;
				pCurrentSatellite->clsActiveDevices[i+1].qQuaternionAttitude._v[2] = Satellite9BQuaternionAttitude_Q2;
				pCurrentSatellite->clsActiveDevices[i+1].qQuaternionAttitude._v[3] = Satellite9BQuaternionAttitude_Q3;

				pCurrentSatellite->clsActiveDevices[i+1].vEulerAttitude._v[0] = Satellite9BEular[0];
				pCurrentSatellite->clsActiveDevices[i+1].vEulerAttitude._v[1] = Satellite9BEular[1];
				pCurrentSatellite->clsActiveDevices[i+1].vEulerAttitude._v[2] = Satellite9BEular[2];

				pCurrentSatellite->clsActiveDevices[i+1].vAngularSpeed._v[0] = Satellite9BAngleSpeed[0];
				pCurrentSatellite->clsActiveDevices[i+1].vAngularSpeed._v[1] = Satellite9BAngleSpeed[1];
				pCurrentSatellite->clsActiveDevices[i+1].vAngularSpeed._v[2] = Satellite9BAngleSpeed[2];

				OrbitUpdate(pCurrentSatellite->clsActiveDevices[i+1].vPosition,pCurrentSatellite->clsActiveDevices[80].pNodeModel->asGeode()->getDrawable(0));
			}		

			if(strcmp(chrDataType,"DecoyAttitudes")==0 && pCurrentManipulator->CurrentChosenDecoyIndex == lCraftNo + 1)
			{
				{
					char chrTemp[100]={'\0'};
					sprintf(chrTemp,"Decoy%d Position:\nRx(m):%lf\nRy(m):%lf\nRz(m):%lf\n",pCurrentManipulator->CurrentChosenDecoyIndex, SatellitePosition_X,SatellitePosition_Y,SatellitePosition_Z);
					WriteToHud(0,0,chrTemp);
				}
				{
					char chrTemp[100]={'\0'};
					sprintf(chrTemp,"Decoy%d Speed:\nVx(m/s):%lf\nVy(m/s):%lf\nVz(m/s):%lf\n",pCurrentManipulator->CurrentChosenDecoyIndex, SatelliteSpeed[0], SatelliteSpeed[1], SatelliteSpeed[2]);
					WriteToHud(0,1,chrTemp);
				}
				{
					char chrTemp[100]={'\0'};
					sprintf(chrTemp,"EularAngle:\nRoll(deg) :%lf\nPitch(deg):%lf\nYew(deg)  :%lf\n",SatelliteEular[0],SatelliteEular[1],SatelliteEular[2]);
					WriteToHud(0,2,chrTemp);
				}
				{	
					char chrTemp[100]={'\0'};
					sprintf(chrTemp,"AngleSpeed:\nX(deg/s):%lf\nY(deg/s):%lf\nZ(deg/s):%lf\n",SatelliteAngleSpeed[0],SatelliteAngleSpeed[1],SatelliteAngleSpeed[2]);
					WriteToHud(0,3,chrTemp);
				}
			}
		}
	}
}


void PickResponse(unsigned char ucDeviceType,unsigned char ucMainBodyIndex,unsigned char ucDeviceIndex)
{
	if (ucDeviceIndex!=0)
	{
		if (pMainBodyInstance[0]->clsActiveDevices[ucDeviceIndex].bHighLight)
		{
			pMainBodyInstance[0]->clsActiveDevices[ucDeviceIndex].CancelHighLight();
		}
		else
		{
			pMainBodyInstance[0]->clsActiveDevices[ucDeviceIndex].bHighLight = true;
		}
	}
	pCurrentManipulator->CurrentChosenDeviceIndex = ucDeviceIndex;

	
}
