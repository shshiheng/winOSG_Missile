#ifndef __CREATE_H
#define __CREATE_H
//#include "Variable.h"
#include <osgViewer/Viewer>
#include <osg/Group>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/Texture2D>
#include <osg/PositionAttitudeTransform>
#include <osg/Matrix>
#include <osgDB/FileUtils>
#include <osg/BlendFunc>
#include <osg/StateSet>
#include <osg/BlendEquation>
#include <osgDB/WriteFile>
#include <osgDB/ReadFile>
#include <osg/State>
#include <osg/Material>
#include <osg/Texture1D>
#include <osg/TexGen>
#include <osg/MatrixTransform>
#include <osg/PolygonStipple>
#include <osg/BlendFunc>
#include <osg/StateSet>
#include <osg/BlendEquation>
//#include "Earth.h"
//#include "GetFormationOrbit.h"
using namespace osg;


osg::Node*  pCreateEarth();
osg::Node* CreatePlusSign();
void pCreateEquator(const char chrMessages[],const unsigned char ucLength);
//osg::Node* pCreateMainOrbit(struct stcOrbitElement *OE,osg::Vec3 vecPointArray[10][360]);
//osg::Node* pCreateFormationOrbit(struct stcOrbitElement OE0,struct stcOrbitElement OE1,osg::Vec4 vColorPlane,osg::Vec4 vColorOrbit,osg::Vec3 vecPointArray[2][360]);
osg::Node* CreateBaseline(osg::Vec3 StartPoint,osg::Vec3 EndPoint,osg::Geometry *pGeometry);
osg::Node* CreateSubstractSign();
osg::Node* CreatePlusSign();
osg::Group* pCreateLight();
osg::Node* pCreateMainOrbit();
osg::StateSet* pCreate1DTextureStateToDecorate(osg::Node* pLoadedModel);
osg::Matrix CalculateRotate(osg::Vec3 SelfAxis,osg::Vec3 vStartPoint,osg::Vec3 vEndPoint);
osg::Matrix CalculateScale(osg::Vec3 vStartPoint,osg::Vec3 vEndPoint);
void OrbitUpdate(osg::Vec3 InVertices,osg::Drawable* drawable);
void RelativeOrbitUpdate(osg::Vec3 InVertices,osg::Drawable* drawable);
osg::Node* pCreateDynamicOrbit();
void CreateProjectionMatrix(osg::Camera *camera);
#endif