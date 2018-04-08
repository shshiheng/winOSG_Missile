//File DigitalGloveUDP contains the definitions of DigitalGloveUDP functions
//File created @ 2009-10-27

#ifndef __PLANET_H
#define __PLANET_H

#include <osg/Node>


class Planet :public osg::NodeCallback
{
private:
	osg::Node* pNodeModel;
	osg::Vec3 vPosition;
	osg::Vec3 vAttitude;
public:
	osg::MatrixTransform* pMatrixTransform;
	bool bVisible;
public:
	Planet();
	osg::Node* CreateModel(const string strPictureName,const double Radius);
	void operator()(const osg::Node *pNode,osg::NodeVisitor *pNodeVisitor);
};

struct StarInfo
{
	double Lattitude;
	double Longitude;
	signed char scLevel;
	unsigned long ulIndex;
};
osg::Node* CreateSkyBackGround(void);
#endif