//File tempManipulator contains the definitions of tempManipulator functions
//File created @ 2009-10-29

#ifndef __TEMPMANIPULATOR_H
#define __TEMPMANIPULATOR_H

#include <iostream>

#include <osg/Transform>
#include <osg/MatrixTransform>
#include <osg/Matrix>
#include <osg/NodeCallback>
#include <osg/Node>
#include <osgGA/CameraManipulator>

#include "Satellite.h"

#define FRAME_TIME_RADIO 0.0167

enum ViewMode{FREE,CIRCULE,FOLLOW_9A,FOLLOW_9B,FOCUS,Watcher_9A};

class SceneManipulator :public osgGA::CameraManipulator
{
public:
	osg::Vec3d vViewPointPosition;
	osg::Vec3d vViewPointPosition1;
	osg::Vec3d vViewDirection;
	osg::Matrixd mManipulator;
	unsigned char CurrentChosenDeviceIndex;
	unsigned char CurrentChosenDecoyIndex;
	ViewMode eViewMode;

	osg::ref_ptr<osg::MatrixTransform> pCrossMatrix;
	osg::Vec3d vCrossPosition;
private:
	osg::Group * pFixedViewObject;
	osgViewer::Viewer *pCurrentViewer;
	osg::Group *pCurrentRoot;
	MainBody *pCenterObject;
	osg::MatrixTransform *pTargetMatrix;
	osg::Quat qRotateQuaternion;
	bool bRotateByQuaternion;
	osg::Vec3d vFocusPoint;
public:
	double CameraLatitude;
	double CameraLongitude;
	double CameraDistence;
private:
	double CameraLatitudeTravelSpeed;
	double CameraLongitudeTravelSpeed;
	double CameraDistenceTravelSpeed;

	bool bViewObject;
	string strManipulatorDiscription;
private:

	inline void CalculateCirculeModePosition();
	inline void CalculateFreeModePosition();
	inline void CalculateFollow9AModePosition();
	inline void CalculateFollow9BModePosition();
	inline void CalculateWatcher9AModePosition();
	inline void CalculateFocusModePosition();
	inline void CameraStep();
	virtual void setByMatrix(const osg::Matrixd& matrix);
	virtual void setByInverseMatrix(const osg::Matrixd& matrix); 
	virtual osg::Matrixd getMatrix(void) const;					
	virtual osg::Matrixd getInverseMatrix(void)const;
public:
	SceneManipulator(osgViewer::Viewer *pViewer,osg::Group *pRoot);
	~SceneManipulator(){};

	void CreateCross();
	virtual bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& us);

	void SetViewObject(osg::MatrixTransform *pTarget);//随后的处理中，应该完成视点的绑定过程
	void SetFollowObject(MainBody *pViewObject);

	void SetViewAreaScale(const double Ratio);
	void SetCameraSpeed(const double LatitudeTravelSpeed,const double LongitudeTravelSpeed,const double DistenceTravelSpeed);

	void IncreasePosition(const osg::Vec3d& vDelta);
	void IncreaseDirection(const osg::Vec3d& vDelta);

	void BackHomeView(void);

	string GetManipulatorDiscription(void);
};

void ManipulatorUpdate(const char chrMessages[],const unsigned char ucLength);

#endif

#ifdef CONNECT_HCI
unsigned short InteractionUseVoiceCommand(const char chrMessages[],const unsigned short ucLength);
unsigned short InteractionUseGestureCommand(const char chrMessages[],const unsigned short ucLength);
#endif
