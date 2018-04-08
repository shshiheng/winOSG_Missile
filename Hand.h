//File Hand contains the definitions of Hand functions
//File created @ 2009-10-29

#ifndef __HAND_H
#define __HAND_H

#pragma once 
#include <osg/NodeCallback>
#include <osg/MatrixTransform>
#include <osg/NodeVisitor>
#include <osg/Node>
#include <osgViewer/Viewer>
#include <osg/Group>
#include <osg/PositionAttitudeTransform>

#include "Manipulator.h"


class HandVisitor :public osg::NodeVisitor
{
public:
	osg::Node* pNodeThumb;
	osg::Node* pNodeThumbFlex;
	osg::Node* pNodeIndex;
	osg::Node* pNodeIndexFlex;
	osg::Node* pNodeMiddle;
	osg::Node* pNodeMiddleFlex;
	osg::Node* pNodeRing;
	osg::Node* pNodeRingFlex;
	osg::Node* pNodeLittle;
	osg::Node* pNodeLittleFlex;
public:
	HandVisitor():osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN){};
	virtual void apply(osg::Node &Node);
};

class Hand :public osg::NodeCallback
{
public:
	osg::Node * pNodeHand;
	osgViewer::Viewer* pCurrentViewer;
	HandVisitor clsHandNodeVisitor;
	bool bFingleState[5];
	osg::Vec3 HandPosition;
	osg::Vec3 HandAttitude;
	float fGloveData[11];
public:
	Hand();
	void ReadSubNode(){pNodeHand->accept(clsHandNodeVisitor);}
	virtual void operator()(osg::Node* pNode,osg::NodeVisitor* pNodeVisitor);
};
class cHandFollowViewpoint : public osg::NodeCallback
{
public:	
	osg::Vec3 viewPosition;
	SceneManipulator *pCurrentManipulator;	
	osgViewer::Viewer* pCurrentViewer;
public:
	cHandFollowViewpoint(osgViewer::Viewer* pViewer,SceneManipulator *pManipulator)
	{
		pCurrentViewer = pViewer;
		pCurrentManipulator = pManipulator;
	};
public:
	virtual void operator()(osg::Node* pNode,osg::NodeVisitor* pNodeVisitor);
};

#endif
