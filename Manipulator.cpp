//File Manipulator contains the realization of Manipulator functions
//File created @ 2009-11-6

#include <iostream>
using namespace std;

#include <osgDB/ReadFile>
#include <osg/CameraNode>
#include <math.h>

#include "Manipulator.h"
#include "Hud.h"
#include "Satellite.h"
#include "CreateDevice.h"
#include "UserMath.h"

#define MAX_MAINBODY_NUMBER 5
#define qInitialQuat osg::Quat(-0.5,-0.5,0.5,0.5)
osg::Matrixd mInveseMatrix;
extern osg::ref_ptr<osgViewer::Viewer> pViewer=new osgViewer::Viewer;
extern SceneManipulator *pCurrentManipulator;
extern MainBody *pCurrentSatellite;
extern MainBody *pMainBodyInstance[MAX_MAINBODY_NUMBER];
osg::ref_ptr<osg::Group> pLight=new osg::Group;
char LightNum = 0;
//Construct Function,Initialize all referenced variables

signed char CheckFireStatus();
SceneManipulator::SceneManipulator(osgViewer::Viewer *pViewer,osg::Group *pRoot)
{
	this->eViewMode = CIRCULE;       //FOCUS   CIRCULE
	this->CameraLongitude = 0.0f;
	this->CameraLatitude = -2.1f;
	this->CameraDistence = 6.6E+8;
	this->vViewPointPosition = osg::Vec3d(0.0, 0.0, 0.0);				//初始位置（视点）
	this->vViewPointPosition1 = osg::Vec3d(0.0, 0.0, 0.0);				//初始位置（视点）
	this->vViewDirection = osg::Vec3d(0.0, 0.0, 0.0);					//初始角度
	this->CameraLatitudeTravelSpeed = 0.0;
	this->CameraLongitudeTravelSpeed = 0.0;
	this->CameraDistenceTravelSpeed = 0.0;
	this->bViewObject = false;
	this->strManipulatorDiscription = "";
	this->mManipulator = osg::Matrix::identity();
	this->pFixedViewObject = new osg::Group ;
	this->pCurrentViewer = pViewer;
	this->pCurrentRoot = pRoot;
	this->pCrossMatrix = new osg::MatrixTransform;
	this->vCrossPosition = osg::Vec3d(0.0,0.0,0.0);
	this->CreateCross();
	this->pCenterObject = NULL;
	this->pTargetMatrix = new osg::MatrixTransform;
	this->CurrentChosenDeviceIndex = 255;
	this->CurrentChosenDecoyIndex = 1;
	this->bRotateByQuaternion = true;
	this->qRotateQuaternion = osg::Quat(0,0,0,1);
}

//functions for Viewport control Matrix
void SceneManipulator::setByMatrix(const osg::Matrixd &matrix)
{
}

void SceneManipulator::setByInverseMatrix(const osg::Matrixd &matrix)
{
}
osg::Matrixd SceneManipulator::getMatrix() const
{
	osg::Matrixd mMatrix;
	osg::Matrixd mOrbit;
	mMatrix.makeRotate(this->vViewDirection._v[0], osg::Vec3d(1.0, 0.0, 0.0), this->vViewDirection._v[1], osg::Vec3d(0.0, 1.0, 0.0), this->vViewDirection._v[2], osg::Vec3d(0.0, 0.0, 1.0));  
	if(this->eViewMode == FOLLOW_9A)
	{	
		double mTranMatrix[3][3] = {0.0};
		double Q[4] = {0.0};
		double R[3] = {pCurrentSatellite->clsActiveDevices[0].vPosition.x(),pCurrentSatellite->clsActiveDevices[0].vPosition.y(),pCurrentSatellite->clsActiveDevices[0].vPosition.z()};
		double V[3] = {pCurrentSatellite->clsActiveDevices[0].vSpeed.x(),pCurrentSatellite->clsActiveDevices[0].vSpeed.y(),pCurrentSatellite->clsActiveDevices[0].vSpeed.z()};
		if(fabs(R[0])<1000)
			;
		else 
		{
			GetRoi(R,V,mTranMatrix);
			QuaternionExtract(mTranMatrix,Q);
			mOrbit.makeRotate(osg::Quat(Q[1],Q[2],Q[3],Q[0]));
		}
		return mMatrix * osg::Matrix::translate(this->vViewPointPosition1) * mOrbit * osg::Matrixd::translate(this->vViewPointPosition);
	}
	else if(this->eViewMode == Watcher_9A)
	{	
		double mTranMatrix[3][3] = {0.0};
		double Q[4] = {0.0};
		double R[3] = {pCurrentSatellite->clsActiveDevices[0].vPosition.x(),pCurrentSatellite->clsActiveDevices[0].vPosition.y(),pCurrentSatellite->clsActiveDevices[0].vPosition.z()};
		double V[3] = {pCurrentSatellite->clsActiveDevices[0].vSpeed.x(),pCurrentSatellite->clsActiveDevices[0].vSpeed.y(),pCurrentSatellite->clsActiveDevices[0].vSpeed.z()};
		if(fabs(R[0])<1000)
			;
		else 
		{
			GetRoi(R,V,mTranMatrix);
			QuaternionExtract(mTranMatrix,Q);
			mOrbit.makeRotate(osg::Quat(Q[1],Q[2],Q[3],Q[0]));
		}
		return mMatrix * osg::Matrix::translate(this->vViewPointPosition1) * mOrbit * osg::Matrixd::translate(this->vViewPointPosition);
	}
	else if(this->eViewMode == FOLLOW_9B)
	{	
		double mTranMatrix[3][3] = {0.0};
		double Q[4] = {0.0};
		double R[3] = {pCurrentSatellite->clsActiveDevices[1].vPosition.x(),pCurrentSatellite->clsActiveDevices[1].vPosition.y(),pCurrentSatellite->clsActiveDevices[1].vPosition.z()};
		double V[3] = {pCurrentSatellite->clsActiveDevices[1].vSpeed.x(),pCurrentSatellite->clsActiveDevices[1].vSpeed.y(),pCurrentSatellite->clsActiveDevices[1].vSpeed.z()};
		if(fabs(R[0])<1000)
			;
		else 
		{
			GetRoi(R,V,mTranMatrix);
			QuaternionExtract(mTranMatrix,Q);
			mOrbit.makeRotate(osg::Quat(Q[1],Q[2],Q[3],Q[0]));
		}
		return mMatrix * osg::Matrix::translate(this->vViewPointPosition1) * mOrbit * osg::Matrixd::translate(this->vViewPointPosition);
	}
	else
		return  mMatrix * osg::Matrixd::translate(this->vViewPointPosition);
}

osg::Matrixd SceneManipulator::getInverseMatrix() const
{
	osg::Matrixd mMatrix;
	osg::Matrixd mOrbit;
	if(this->eViewMode == FOCUS)
	{
		/*double xLen = this->vFocusPoint._v[0]-this->vViewPointPosition._v[0];
		double yLen = this->vFocusPoint._v[1]-this->vViewPointPosition._v[1];
		double zLen = this->vFocusPoint._v[2]-this->vViewPointPosition._v[2];
		double length = sqrt(xLen*xLen+yLen*yLen+zLen*zLen);
		double angleX = xLen/length;
		double angleY = yLen/length;
		double angleZ = zLen/length;*/


		mMatrix = osg::Matrixd::rotate(osg::PI_2,osg::X_AXIS)*osg::Matrixd::rotate(this->vViewDirection._v[0],osg::X_AXIS)*osg::Matrixd::rotate(this->vViewDirection._v[2],osg::Z_AXIS);
	}
	else
		mMatrix.makeRotate(this->vViewDirection._v[0], osg::Vec3d(1.0, 0.0, 0.0), this->vViewDirection._v[1], osg::Vec3d(0.0, 1.0, 0.0), this->vViewDirection._v[2], osg::Vec3d(0.0, 0.0, 1.0)); 

	if(this->eViewMode == FOLLOW_9A)
	{
		double mTranMatrix[3][3] = {0.0};
		double Q[4] = {0.0};
		double R[3] = {pCurrentSatellite->clsActiveDevices[0].vPosition.x(),pCurrentSatellite->clsActiveDevices[0].vPosition.y(),pCurrentSatellite->clsActiveDevices[0].vPosition.z()};
		double V[3] = {pCurrentSatellite->clsActiveDevices[0].vSpeed.x(),pCurrentSatellite->clsActiveDevices[0].vSpeed.y(),pCurrentSatellite->clsActiveDevices[0].vSpeed.z()};
		if(fabs(R[0])<10)
			;
		else 
		{
			GetRoi(R,V,mTranMatrix);
			QuaternionExtract(mTranMatrix,Q);
			mOrbit.makeRotate(osg::Quat(Q[1],Q[2],Q[3],Q[0]));
		}
		return mInveseMatrix = osg::Matrixd::inverse(mMatrix * osg::Matrix::translate(this->vViewPointPosition1) * mOrbit * osg::Matrixd::translate(this->vViewPointPosition));
	}
	else if(this->eViewMode == Watcher_9A)
	{
		double mTranMatrix[3][3] = {0.0};
		double Q[4] = {0.0};
		double R[3] = {pCurrentSatellite->clsActiveDevices[0].vPosition.x(),pCurrentSatellite->clsActiveDevices[0].vPosition.y(),pCurrentSatellite->clsActiveDevices[0].vPosition.z()};
		double V[3] = {pCurrentSatellite->clsActiveDevices[0].vSpeed.x(),pCurrentSatellite->clsActiveDevices[0].vSpeed.y(),pCurrentSatellite->clsActiveDevices[0].vSpeed.z()};

		return mInveseMatrix = osg::Matrixd::inverse(mMatrix * osg::Matrix::translate(this->vViewPointPosition1) * mOrbit * osg::Matrixd::translate(this->vViewPointPosition));
	}
	else if(this->eViewMode == FOLLOW_9B)
	{
		double mTranMatrix[3][3] = {0.0};
		double Q[4] = {0.0};
		double R[3] = {pCurrentSatellite->clsActiveDevices[1].vPosition.x(),pCurrentSatellite->clsActiveDevices[1].vPosition.y(),pCurrentSatellite->clsActiveDevices[1].vPosition.z()};
		double V[3] = {pCurrentSatellite->clsActiveDevices[1].vSpeed.x(),pCurrentSatellite->clsActiveDevices[1].vSpeed.y(),pCurrentSatellite->clsActiveDevices[1].vSpeed.z()};
		if(fabs(R[0])<10)
			;
		else 
		{
			GetRoi(R,V,mTranMatrix);
			QuaternionExtract(mTranMatrix,Q);
			mOrbit.makeRotate(osg::Quat(Q[1],Q[2],Q[3],Q[0]));
		}
		return mInveseMatrix =osg::Matrixd::inverse(mMatrix * osg::Matrix::translate(this->vViewPointPosition1) * mOrbit * osg::Matrixd::translate(this->vViewPointPosition));
	}
	else
		return  mInveseMatrix =osg::Matrixd::inverse(mMatrix * osg::Matrixd::translate(this->vViewPointPosition));
}

//set the ceter object to focus
void SceneManipulator::SetViewObject(osg::MatrixTransform *pMatrix)
{
	if(pMatrix != NULL)
	{
#ifdef DEBUG
		cout<<"Set a new Focus Target for viewer "<<endl;
#endif
		this->pTargetMatrix = pMatrix;
		this->eViewMode = FOLLOW_9A;

	}
	else
	{
#ifdef DEBUG
		cout<<"Back to CIRCULE Mode"<<endl;
#endif
		this->eViewMode = CIRCULE;
	}
}

void SceneManipulator::SetFollowObject(MainBody *pTarget)
{
	if(pTarget == NULL)
	{
#ifdef DEBUG
		cout<<"Back to CIRCULE Mode"<<endl;
#endif
		this->eViewMode = CIRCULE;
	}
	else
	{
#ifdef DEBUG
		cout<<"Set a new following Target for viewer "<<endl;
#endif
		this->pCenterObject = pTarget;
		this->eViewMode = FOLLOW_9A;
	}
}

inline void SceneManipulator::CalculateCirculeModePosition()
{//Code for Circule Mode
	this->bRotateByQuaternion = false;
	this->eViewMode = CIRCULE;
	this->vViewPointPosition.x() = this->CameraDistence*cos(this->CameraLatitude)*sin(this->CameraLongitude);
	this->vViewPointPosition.y() = -this->CameraDistence*cos(this->CameraLatitude)*cos(this->CameraLongitude);
	this->vViewPointPosition.z() = this->CameraDistence*sin(this->CameraLatitude);			
	this->vViewDirection._v[0] = osg::PI_2-this->CameraLatitude;
	this->vViewDirection._v[1] = 0;
	this->vViewDirection._v[2] = this->CameraLongitude;
}

inline void SceneManipulator::CalculateFreeModePosition()
{//Code for Free Mode,camera can move freely
	this->vViewDirection._v[0] = this->CameraLatitude;
	this->vViewDirection._v[1] = 0;
	this->vViewDirection._v[2] = this->CameraLongitude;
}

inline void SceneManipulator::CalculateFollow9AModePosition()
{//Code for FOLLOW Mode
	if(this->pCenterObject == NULL)
	{
#ifdef DEBUG
		cout<<"Invaild Following Target, Please Check"<<endl;
#endif
		this->eViewMode = CIRCULE;
	}
	else
	{
		this->eViewMode = FOLLOW_9A;
		this->bRotateByQuaternion = true;
		osg::Vec3d vBasePosition = this->pCenterObject->clsActiveDevices[0].vPosition;
		osg::Vec3d vBaseAttitude = this->pCenterObject->clsActiveDevices[0].vEulerAttitude;
		this->vViewPointPosition1._v[0] = this->CameraDistence*cos(this->CameraLatitude)*sin(this->CameraLongitude);
		this->vViewPointPosition1._v[1] = - this->CameraDistence*cos(this->CameraLatitude)*cos(this->CameraLongitude);
		this->vViewPointPosition1._v[2] = this->CameraDistence*sin(this->CameraLatitude);


		this->vViewPointPosition._v[0] = vBasePosition._v[0];
		this->vViewPointPosition._v[1] = vBasePosition._v[1];
		this->vViewPointPosition._v[2] = vBasePosition._v[2];


		this->vViewDirection._v[0] = osg::PI_2-this->CameraLatitude;
		this->vViewDirection._v[1] = 0;
		this->vViewDirection._v[2] = this->CameraLongitude;

		this->qRotateQuaternion._v[0] = this->pCenterObject->clsActiveDevices[0].qQuaternionAttitude[0];
		this->qRotateQuaternion._v[1] = this->pCenterObject->clsActiveDevices[0].qQuaternionAttitude[1];
		this->qRotateQuaternion._v[2] = this->pCenterObject->clsActiveDevices[0].qQuaternionAttitude[2];
		this->qRotateQuaternion._v[3] = this->pCenterObject->clsActiveDevices[0].qQuaternionAttitude[3];

	}
	this->eViewMode = FOLLOW_9A;

}

inline void SceneManipulator::CalculateWatcher9AModePosition()
{//Code for Watcher9A Mode
	if(this->pCenterObject == NULL)
	{
#ifdef DEBUG
		cout<<"Invaild Following Target, Please Check"<<endl;
#endif
		this->eViewMode = CIRCULE;
	}
	else
	{
		this->eViewMode = Watcher_9A;
		this->bRotateByQuaternion = false;
		osg::Vec3d vBasePosition = this->pCenterObject->clsActiveDevices[0].vPosition;
		osg::Vec3d vBaseAttitude = this->pCenterObject->clsActiveDevices[0].vEulerAttitude;
		this->vViewPointPosition1._v[0] = this->CameraDistence*cos(this->CameraLatitude)*sin(this->CameraLongitude);
		this->vViewPointPosition1._v[1] = - this->CameraDistence*cos(this->CameraLatitude)*cos(this->CameraLongitude);
		this->vViewPointPosition1._v[2] = this->CameraDistence*sin(this->CameraLatitude);


		this->vViewPointPosition._v[0] = vBasePosition._v[0];
		this->vViewPointPosition._v[1] = vBasePosition._v[1];
		this->vViewPointPosition._v[2] = vBasePosition._v[2];


		this->vViewDirection._v[0] = osg::PI_2-this->CameraLatitude;
		this->vViewDirection._v[1] = 0;
		this->vViewDirection._v[2] = this->CameraLongitude;


	}
	this->eViewMode = Watcher_9A;
}

inline void SceneManipulator::CalculateFollow9BModePosition()
{//Code for FOLLOW Mode
	if(this->pCenterObject == NULL)
	{
#ifdef DEBUG
		cout<<"Invaild Following Target, Please Check"<<endl;
#endif
		this->eViewMode = CIRCULE;
	}
	else
	{
		this->eViewMode = FOLLOW_9B;
		this->bRotateByQuaternion = true;
		osg::Vec3d vBasePosition = this->pCenterObject->clsActiveDevices[1].vPosition;
		osg::Vec3d vBaseAttitude = this->pCenterObject->clsActiveDevices[1].vEulerAttitude;
		this->vViewPointPosition1._v[0] = this->CameraDistence*cosf(this->CameraLatitude)*sinf(this->CameraLongitude);
		this->vViewPointPosition1._v[1] = - this->CameraDistence*cosf(this->CameraLatitude)*cosf(this->CameraLongitude);
		this->vViewPointPosition1._v[2] = this->CameraDistence*sinf(this->CameraLatitude);


		this->vViewPointPosition._v[0] = vBasePosition._v[0];
		this->vViewPointPosition._v[1] = vBasePosition._v[1];
		this->vViewPointPosition._v[2] = vBasePosition._v[2];


		this->vViewDirection._v[0] = osg::PI_2-this->CameraLatitude;
		this->vViewDirection._v[1] = 0;
		this->vViewDirection._v[2] = this->CameraLongitude;

		this->qRotateQuaternion._v[0] = this->pCenterObject->clsActiveDevices[1].qQuaternionAttitude[0];
		this->qRotateQuaternion._v[1] = this->pCenterObject->clsActiveDevices[1].qQuaternionAttitude[1];
		this->qRotateQuaternion._v[2] = this->pCenterObject->clsActiveDevices[1].qQuaternionAttitude[2];
		this->qRotateQuaternion._v[3] = this->pCenterObject->clsActiveDevices[1].qQuaternionAttitude[3];

	}
	this->eViewMode = FOLLOW_9B;
}

inline void SceneManipulator::CalculateFocusModePosition()
{//code for calculate focus mode position
	osg::MatrixList mList = this->pTargetMatrix->getWorldMatrices();
	osg::Matrix finalMatrix = osg::Matrix::identity();
	for(vector <osg::Matrix>::iterator it = mList.begin();it != mList.end(); ++it)
	{
		osg::Matrix mat = *it;
		finalMatrix = finalMatrix*mat;
	}
	this->vFocusPoint = finalMatrix.getTrans();

	this->vViewPointPosition._v[0] = this->vFocusPoint._v[0] + this->CameraDistence*cosf(this->CameraLatitude)*sinf(this->CameraLongitude);
	this->vViewPointPosition._v[1] = this->vFocusPoint._v[1] - this->CameraDistence*cosf(this->CameraLatitude)*cosf(this->CameraLongitude);
	this->vViewPointPosition._v[2] = this->vFocusPoint._v[2] + this->CameraDistence*sinf(this->CameraLatitude);

	this->vViewDirection._v[0] = - this->CameraLatitude;
	this->vViewDirection._v[1] = 0;
	this->vViewDirection._v[2] = this->CameraLongitude;
}

inline void SceneManipulator::CameraStep()
{//Code to change Camera Params every frame
	this->CameraLatitude += this->CameraLatitudeTravelSpeed;
	this->CameraLongitude += this->CameraLongitudeTravelSpeed;
	this->CameraDistence += this->CameraDistenceTravelSpeed;



	//if(this->CameraLatitude >= 2*osg::PI)
	//	this->CameraLatitude = this->CameraLatitude - 2*osg::PI;
	//else if(this->CameraLatitude < -2*osg::PI)
	//	this->CameraLatitude = this->CameraLatitude + 2*osg::PI;
	//if(this->CameraLongitude >= 2*osg::PI)
	//	this->CameraLongitude = this->CameraLongitude - 2*osg::PI;
	//else if(this->CameraLongitude < -2*osg::PI)
	//	this->CameraLongitude = this->CameraLongitude + 2*osg::PI;
}

bool SceneManipulator::handle(const osgGA::GUIEventAdapter &ea, osgGA::GUIActionAdapter &us)
{
	static bool bChangeFireEffect = false;
	if(eViewMode != FOCUS)
	{
	}

	switch(ea.getEventType())
	{
	case osgGA::GUIEventAdapter::FRAME:
		{
			this->CameraStep();
			if(this->eViewMode == CIRCULE)		//Code calculate new position of camera
				this->CalculateCirculeModePosition();	
			else if(this->eViewMode == FOLLOW_9A)
				this->CalculateFollow9AModePosition();
			else if(this->eViewMode == FOLLOW_9B)
				this->CalculateFollow9BModePosition();
			else if(this->eViewMode == FREE)
				this->CalculateFreeModePosition();
			else if(this->eViewMode == FOCUS)
				this->CalculateFocusModePosition();
			else if(this->eViewMode == Watcher_9A)
				this->CalculateWatcher9AModePosition();

			this->pCrossMatrix->setMatrix(osg::Matrix::translate(pCurrentSatellite->clsActiveDevices[1].vPosition));
			CheckFireStatus();				
			return false;
		}		
	case osgGA::GUIEventAdapter::LEFT_MOUSE_BUTTON:
		{		
			double MouseXPosition = ea.getX();
			double MouseYPosition = ea.getY();
			Pick(MouseXPosition,MouseYPosition);

			cout<<"xPos:  "<<MouseXPosition<<endl;
			cout<<"yPos:  "<<MouseYPosition<<endl;

			return false;
		}	
	case osgGA::GUIEventAdapter::KEYDOWN:
		{		
			if (ea.getKey() == 'A'||ea.getKey() == 'a')		//Left 
			{	
				this->CameraLongitude += osg::DegreesToRadians(5.0f);
				return false;
			}	
			if(ea.getKey() == 'W'||ea.getKey() == 'w')		//Up
			{	
				this->CameraLatitude += osg::DegreesToRadians(5.0f);	
				return false;
			}	
			if (ea.getKey() == 'D'||ea.getKey() == 'd')		//Right
			{	
				this->CameraLongitude -= osg::DegreesToRadians(5.0f);				
				return false;
			}	
			if(ea.getKey() == 'S'||ea.getKey() == 's')		//Down
			{	
				this->CameraLatitude -= osg::DegreesToRadians(5.0f);
				return false;
			}	

			if(ea.getKey() == 'Q'||ea.getKey() == 'q')	//Q
			{	
				if(this->CameraDistence < 1000000)
					this->CameraDistence += 10000;

				else if(this->CameraDistence >= 1000000 && this->CameraDistence < 10000000)
					this->CameraDistence += 500000;
				else
					this->CameraDistence += 50000000;
				return false;			
			}	
			if (ea.getKey() == 'E'||ea.getKey() == 'e')	//E
			{	
				if(this->CameraDistence < 1000000)
					this->CameraDistence -= 10000;

				else if(this->CameraDistence >= 1000000 && this->CameraDistence <= 10000000)
					this->CameraDistence -= 500000;
				else if(this->CameraDistence > 10000000 && this->CameraDistence <= 100000000)
					this->CameraDistence -= 5000000;
				else
					this->CameraDistence -= 50000000;
				return false;			
			}	
			if (ea.getKey() == 'N'||ea.getKey() == 'n')			
			{	
				this->bRotateByQuaternion = false;
				this->eViewMode = CIRCULE;
				this->CameraLatitude = osg::DegreesToRadians(31.5f);
				this->CameraLongitude = osg::DegreesToRadians(15.5f);
				this->CameraDistence = 6.6E+8;
				return false;			
			}	
			if(ea.getKey() == 'M'||ea.getKey() == 'm')
			{	
				this->eViewMode  = FOLLOW_9A;
				this->CameraLongitude = 0.0f;
				this->CameraLatitude = -2.1f;
				this->CameraDistence = 6.6E+8;
				return false;
			}	
			if(ea.getKey() == 'B'||ea.getKey() == 'b')
			{	
				this->eViewMode  = FOLLOW_9B;
				this->CameraLongitude = 0.0f;
				this->CameraLatitude = -2.1f;
				this->CameraDistence = 6.6E+8;
				return false;
			}	
			if(ea.getKey() == '1')
			{	
				this->CurrentChosenDecoyIndex = 1;
				return false;
			}
			if(ea.getKey() == '2')
			{	
				this->CurrentChosenDecoyIndex = 2;
				return false;
			}
			if(ea.getKey() == '3')
			{	
				this->CurrentChosenDecoyIndex = 3;
				return false;
			}
			if(ea.getKey() == '4')
			{	
				this->CurrentChosenDecoyIndex = 4;
				return false;
			}
			if(ea.getKey() == '5')
			{	
				this->CurrentChosenDecoyIndex = 5;
				return false;
			}
			if(ea.getKey() == '6')
			{	
				this->CurrentChosenDecoyIndex = 6;
				return false;
			}
			if(ea.getKey() == '7')
			{	
				this->CurrentChosenDecoyIndex = 7;
				return false;
			}
			if(ea.getKey() == '8')
			{	
				this->CurrentChosenDecoyIndex = 8;
				return false;
			}
			if(ea.getKey() == '9')
			{	
				this->CurrentChosenDecoyIndex = 9;
				return false;
			}
			if(ea.getKey() == '0')
			{	
				this->CurrentChosenDecoyIndex = 10;
				return false;
			}
			if(ea.getKey()=='L'||ea.getKey()=='l')
			{	
				if(LightNum == 0)
				{
					pLight = pCreateLight();
					LightNum++;
					this->pCurrentRoot->addChild(pLight.get()); 
				}

			}
			if(ea.getKey()=='K'||ea.getKey()=='k')
			{	
				if(LightNum == 1)
				{
					pLight->removeChild(pLight->getChild(0));
					LightNum = 0;
				}
			}
			if(ea.getKey()==0xFF1B)
			{	
				exit(0);
			}

			cout<<"Cam Latitude:"<<this->CameraLatitude<<"\tCam Longitude"<<this->CameraLongitude<<"\tCam Disitence"<<this->CameraDistence<<endl;

			return false;
		}		
	default :	
		return false;
	}			
}				


void SceneManipulator::IncreasePosition(const osg::Vec3d &vDelta)
{				
	this->vViewPointPosition += vDelta;
}				

void SceneManipulator::IncreaseDirection(const osg::Vec3d &vDelta)
{				
	this->vViewDirection += vDelta;
}				

void SceneManipulator::SetViewAreaScale(const double Ratio)
{				
}				

void SceneManipulator::SetCameraSpeed(const double LatitudeTravelSpeed, const double LongitudeTravelSpeed, const double DistenceTravelSpeed)
{				
	this->CameraLatitudeTravelSpeed = LatitudeTravelSpeed;
	this->CameraLongitudeTravelSpeed = LongitudeTravelSpeed;
	this->CameraDistenceTravelSpeed = DistenceTravelSpeed;
}				

void SceneManipulator::BackHomeView()
{				
}				

string SceneManipulator::GetManipulatorDiscription()
{				
	return "";	
}				

void SceneManipulator::CreateCross()
{				
	osg::Camera* camera = new osg::CameraNode;

	camera->setProjectionMatrix(osg::Matrix::ortho2D(0,1024,0,768));	//设置透视矩阵
	camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);	
	camera->setViewMatrix(osg::Matrix::identity());						//得到默认设置	
	camera->setClearMask(GL_DEPTH_BUFFER_BIT);							//设置背景为透明，否则的话可以设置ClearColor   
	camera->setRenderOrder(osg::CameraNode::POST_RENDER);				//设置渲染顺序，必须在最后渲染

	osg::ref_ptr<osg::Node> pCrossNode = osgDB::readNodeFile("DirectionIndicator.3DS");
	//this->pCrossMatrix->setMatrix(osg::Matrix::scale(1000000,1000000,1000000));
	this->pCrossMatrix->addChild(pCrossNode.get());
	this->pCrossMatrix->setNodeMask(0);
	camera->addChild(this->pCrossMatrix.get());

	this->pCurrentRoot->addChild(camera);

}

unsigned short InteractionUseVoiceCommand(const char chrMessages[],const unsigned short ucLength)
{

	char chrManipulatorMessages[4] = {0};
	memcpy(chrManipulatorMessages,chrMessages,ucLength);

	if(1==1)
	{
		switch(chrManipulatorMessages[0])
		{
		case 1:
			pCurrentManipulator->SetCameraSpeed(0.0,osg::DegreesToRadians(0.1f),0.0);
			cout<<"正在向左"<<endl;
			break;
		case 2:
			pCurrentManipulator->SetCameraSpeed(0.0,-osg::DegreesToRadians(0.1f),0.0);
			cout<<"正在向右"<<endl;
			break;
		case 3:
			pCurrentManipulator->SetCameraSpeed(-osg::DegreesToRadians(0.1f),0.0,0.0);
			cout<<"正在向上"<<endl;
			break;
		case 4:
			pCurrentManipulator->SetCameraSpeed(osg::DegreesToRadians(0.1f),0.0,0.0);
			cout<<"正在向下"<<endl;
			break;
		case 5:
			pCurrentManipulator->SetCameraSpeed(0.0,0.0,800.0);
			cout<<"正在缩小"<<endl;
			break;
		case 6:
			pCurrentManipulator->SetCameraSpeed(0.0,0.0,-800.0);
			cout<<"正在放大"<<endl;
			break;
		case 7:
			pCurrentManipulator->SetCameraSpeed(0.0,0.0,0.0);
			cout<<"处于静止状态"<<endl;
			break;
		case 9:
			cout<<"你好！欢迎来到人机交互三维仿真世界，请问有什么可以帮您的？"<<endl;
			break;
		default :
			break;
		}

	}
	return 0;
}

unsigned short InteractionUseGestureCommand(const char chrMessages[],const unsigned short ucLength)
{
	//激活的功能
	char chrGestureMessages[4] = {0};
	memcpy(chrGestureMessages,chrMessages,ucLength);

	if(1==1)
	{
		switch(chrGestureMessages[0])
		{
		case 1:
			pCurrentManipulator->SetCameraSpeed(0.0,osg::DegreesToRadians(0.1f),0.0);
			cout<<"正在向左"<<endl;
			break;
		case 2:
			pCurrentManipulator->SetCameraSpeed(0.0,-osg::DegreesToRadians(0.1f),0.0);
			cout<<"正在向右"<<endl;
			break;
		case 3:
			pCurrentManipulator->SetCameraSpeed(-osg::DegreesToRadians(0.1f),0.0,0.0);
			cout<<"正在向上"<<endl;
			break;
		case 4:
			pCurrentManipulator->SetCameraSpeed(osg::DegreesToRadians(0.1f),0.0,0.0);
			cout<<"正在向下"<<endl;
			break;
		case 5:
			pCurrentManipulator->SetCameraSpeed(0.0,0.0,800.0);
			cout<<"正在缩小"<<endl;
			break;
		case 6:
			pCurrentManipulator->SetCameraSpeed(0.0,0.0,-800.0);
			cout<<"正在放大"<<endl;
			break;
		case 7:
			pCurrentManipulator->SetCameraSpeed(0.0,0.0,0.0);
			cout<<"处于静止状态"<<endl;
			break;
		case 9:
			cout<<"你好！欢迎来到人机交互三维仿真世界，请问有什么可以帮您的？"<<endl;
			break;
		default :
			break;
		}

	}
	return 0;
}