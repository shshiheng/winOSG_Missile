#pragma once

#include <vector>
#include <string>

#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Group>
#include <osg/MatrixTransform>
#include <osg/BlendFunc>
#include <osgViewer/Viewer>
#include <osgGA/GUIEventAdapter>
#include <osgGA/GUIActionAdapter>
#include <osgFX/Effect>
#include <osg/Depth>
#include <osg/CameraNode>
#include <osgText/Text>
#include <osgDB/ReadFile>
#include <osg/StateSet>
#include <osg/Node>
#include <osg/ShapeDrawable>

#include "Hud.h"


#define MAX_REGION_COUNT 5
#define MAX_SCROLL_LINE 6

struct Region
{
	osgText::Text *pShowText;
	unsigned char ucMaxCharacters;
	string strContent;
	bool bUpdated;
};

class HUD
{
public:
	Region stcOutputRegion[MAX_REGION_COUNT];
	osg::Group *pGroupHud;	
public:
	bool bVisible;
	bool bAutoWordWrap;
public:
	HUD(){
		for(int iCount=0;iCount<MAX_REGION_COUNT;iCount++)
		{
			this->stcOutputRegion[iCount].pShowText = new osgText::Text;
			this->stcOutputRegion[iCount].pShowText->setDataVariance(osg::Object::DYNAMIC);
			this->stcOutputRegion[iCount].ucMaxCharacters = 1000;
			this->stcOutputRegion[iCount].strContent = "";
			this->stcOutputRegion[iCount].bUpdated = false;
		}
		//set Text Size,font,color,region position,static text contents
	};
	void CreateHud(osg::Group* pRoot,const osg::Vec3 vStartPoint=osg::Vec3(0,768,0),const osg::Vec3 vEndPoint = osg::Vec3(200,668,0));
};

class ScrollHud
{
public:
	osgText::Text *pTextScrollHud[MAX_SCROLL_LINE];
	osg::Geode *pGeodeScrollHud;
	string strShowContent[MAX_SCROLL_LINE];
	bool bUpdated;
public:
	bool bVisible;
public:
	ScrollHud();
	void CreateScrollHud(osg::Group* pRoot,const osg::Vec3 vStartPoint=osg::Vec3(0,667,0),const osg::Vec3 vEndPoint = osg::Vec3(400,467,0));
	void AddString(string strShowContent);
};

class HUDUpdate :public osg::NodeCallback
{
public:
	virtual void operator()(osg::Node* node,osg::NodeVisitor* nv);

};

osg::Geometry* CreateShadow(const osg::Vec3 vStartPoint,const osg::Vec3 vEndPoint)
{
	osg::Geometry *pGeometry =new osg::Geometry;
	osg::Vec3Array  *pVertexes = new osg::Vec3Array;

	pVertexes->push_back(vStartPoint);
	pVertexes->push_back(osg::Vec3(vStartPoint._v[0],vEndPoint._v[1],0));
    pVertexes->push_back(vEndPoint);
	pVertexes->push_back(osg::Vec3(vEndPoint._v[0],vStartPoint._v[1],0));

	pGeometry->setVertexArray(pVertexes);
    osg::Vec3Array* normals = new osg::Vec3Array;
    normals->push_back(osg::Vec3(0.0f,0.0f,1.0f));
    pGeometry->setNormalArray(normals);
    pGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);

	osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(0,0,0,0.0));
    pGeometry->setColorArray(colors);
    pGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);

    pGeometry->addPrimitiveSet(new osg::DrawArrays(GL_QUADS,0,4));

    osg::StateSet* stateset = pGeometry->getOrCreateStateSet();
    stateset->setMode(GL_BLEND,osg::StateAttribute::ON);
	stateset->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);

    return pGeometry;
}

void HUD::CreateHud(osg::Group* pRoot,const osg::Vec3 vStartPoint,const osg::Vec3 vEndPoint)
{
	osg::Geode* pGeode = new osg::Geode();

	osg::StateSet* pStateset=pGeode->getOrCreateStateSet();
	pStateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);	//close Light
	pStateset->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);	// close depth

	pGeode->addDrawable(CreateShadow(vStartPoint,vEndPoint));

	std::string Font("fonts/SIMSUN.TTC");
	osg::Vec4 vColor = osg::Vec4(0,1,1,1);

	float fCharacterSize = 15.0f;
	//float fDistanceToBound = 10;
								
	for(int iIndex = 0;iIndex < MAX_REGION_COUNT;iIndex++)
	{
		pGeode->addDrawable(stcOutputRegion[iIndex].pShowText);
		this->stcOutputRegion[iIndex].pShowText->setMaximumWidth(stcOutputRegion[iIndex].ucMaxCharacters);
		this->stcOutputRegion[iIndex].pShowText->setCharacterSize(fCharacterSize);
		this->stcOutputRegion[iIndex].pShowText->setFont(Font);
		this->stcOutputRegion[iIndex].pShowText->setColor(vColor);
	}

	{
		osg::Camera* camera = new osg::CameraNode;
		camera->setProjectionMatrix(osg::Matrix::ortho2D(0,1024,0,768));	//设置透视矩阵
		camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);	
		camera->setViewMatrix(osg::Matrix::identity());						//得到默认设置	
		camera->setRenderOrder(osg::CameraNode::POST_RENDER);				//设置渲染顺序，必须在最后渲染
		camera->addChild(pGeode);
		camera->setClearMask(GL_DEPTH_BUFFER_BIT);							//设置背景为透明，否则的话可以设置ClearColor   
		pRoot->addChild(camera);
	}
	//return camera;
};

ScrollHud::ScrollHud()
{
	this->pGeodeScrollHud=new osg::Geode;
	for(unsigned char ucIndex=0;ucIndex<MAX_SCROLL_LINE;ucIndex++)
	{
		this->pTextScrollHud[ucIndex]=new osgText::Text;
		this->pTextScrollHud[ucIndex]->setDataVariance(osg::Object::DYNAMIC);
		this-> pTextScrollHud[ucIndex]->setFont("fonts/SIMSUN.TTC");
		this->pTextScrollHud[ucIndex]->setCharacterSize(30.0f);
		this->strShowContent[ucIndex]="NULL";
	}
	bUpdated = false;
}

void ScrollHud::CreateScrollHud(osg::Group* pRoot,const osg::Vec3 vStartPoint,const osg::Vec3 vEndPoint)
{
	double hight=vEndPoint._v[1]-vStartPoint._v[1];
	osg::Vec3 delta=osg::Vec3(0,hight/MAX_SCROLL_LINE,0);
	osg::Vec3 vPosition = vStartPoint;
	float fDistanceToBound = 10;
	vPosition._v[0] = vPosition._v[0] + fDistanceToBound;
	vPosition._v[1] = vPosition._v[1] - fDistanceToBound -20;

	osg::Camera* camera = new osg::CameraNode;
	osg::Geode * pGeode=new osg::Geode;

	osg::StateSet* pStateset=pGeode->getOrCreateStateSet();
	pStateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);	//close Light
	pStateset->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);	// close depth

	pGeode->addDrawable(CreateShadow(vStartPoint,vEndPoint));
    
	camera->setProjectionMatrix(osg::Matrix::ortho2D(0,1024,0,768));	//设置透视矩阵   
	camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);	
    camera->setViewMatrix(osg::Matrix::identity());						//得到默认设置	
    camera->setClearMask(GL_DEPTH_BUFFER_BIT);							//设置背景为透明，否则的话可以设置ClearColor    
	camera->setRenderOrder(osg::CameraNode::POST_RENDER);				//设置渲染顺序，必须在最后渲染
	camera->addChild(pGeode);

	for(unsigned char ucIndex=0;ucIndex<MAX_SCROLL_LINE;ucIndex++)
	{
		this->pTextScrollHud[ucIndex]->setPosition(vPosition);
		vPosition = vPosition + delta;
		pGeode->addDrawable(this->pTextScrollHud[ucIndex]);
	}
	pRoot->addChild(camera);
}




#define MAX_HUD_NUMBER 5
#define MAX_SCROLLHUD_NUMBER 3

HUD clsHudInstance[MAX_HUD_NUMBER];
ScrollHud clsScrollHudInstance[MAX_SCROLLHUD_NUMBER];

void HUDUpdate::operator ()(osg::Node *node, osg::NodeVisitor *nv)
{
	for(unsigned char ucHUDIndex = 0;ucHUDIndex < MAX_HUD_NUMBER;ucHUDIndex ++)
	{
		for(unsigned char ucRegionIndex = 0;ucRegionIndex < MAX_REGION_COUNT;ucRegionIndex ++)
		{
			if(clsHudInstance[ucHUDIndex].stcOutputRegion[ucRegionIndex].bUpdated == true)
			{
				clsHudInstance[ucHUDIndex].stcOutputRegion[ucRegionIndex].pShowText->setText(clsHudInstance[ucHUDIndex].stcOutputRegion[ucRegionIndex].strContent);
				clsHudInstance[ucHUDIndex].stcOutputRegion[ucRegionIndex].bUpdated = false;
			}
		}
	}

	for(unsigned char ucScrollHUDIndex = 0;ucScrollHUDIndex < MAX_SCROLLHUD_NUMBER; ucScrollHUDIndex++)
	{
		if(clsScrollHudInstance[ucScrollHUDIndex].bUpdated == true)
		{
			for(unsigned char ucLineIndex = 0;ucLineIndex < MAX_SCROLL_LINE;ucLineIndex ++)
			{
				clsScrollHudInstance[ucScrollHUDIndex].pTextScrollHud[ucLineIndex]->setText(clsScrollHudInstance[ucScrollHUDIndex].strShowContent[ucLineIndex]);
			}
			clsScrollHudInstance[ucScrollHUDIndex].bUpdated = false;
		}
	}
}

void WriteToHud(const unsigned char ucHudIndex,const unsigned char ucRegionIndex,const string strOutput)
{
	clsHudInstance[ucHudIndex].stcOutputRegion[ucRegionIndex].bUpdated = true;
	clsHudInstance[ucHudIndex].stcOutputRegion[ucRegionIndex].strContent = strOutput;
}

void WriteToScrollHud(const unsigned char ucScrollHudIndex,string strShowContent)
{
	clsScrollHudInstance[ucScrollHudIndex].bUpdated = true;
	for(unsigned char ucLineIndex = 0;ucLineIndex < MAX_SCROLL_LINE;ucLineIndex ++)
	{
		clsScrollHudInstance[ucScrollHudIndex].strShowContent[ucLineIndex] = clsScrollHudInstance[ucScrollHudIndex].strShowContent[ucLineIndex + 1];
	}
	clsScrollHudInstance[ucScrollHudIndex].strShowContent[MAX_SCROLL_LINE - 1] = strShowContent;
}


/*****************************************************************

此函数允许程序设计人员安排HUD显示：
	设计内容：HUD及ScrollHUD位置，字体，大小；


*****************************************************************/

void AddHUDToScene(osg::Group* pRoot)
{
	osg::ref_ptr<osg::Group> pHudGroup = new osg::Group;

	{
		clsHudInstance[0].CreateHud(pHudGroup.get());

		osg::Vec3 vPosition(10,745,0);				//  code here change position of text
		clsHudInstance[0].stcOutputRegion[0].pShowText->setPosition(vPosition);	
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[0].stcOutputRegion[1].pShowText->setPosition(vPosition);
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[0].stcOutputRegion[2].pShowText->setPosition(vPosition);
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[0].stcOutputRegion[3].pShowText->setPosition(vPosition);
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[0].stcOutputRegion[4].pShowText->setPosition(vPosition);
	}

	{
		clsHudInstance[1].CreateHud(pHudGroup.get(),osg::Vec3(824,768,0),osg::Vec3(1024,668,0));

		osg::Vec3 vPosition(834,745,0);				//  code here change position of text
		clsHudInstance[1].stcOutputRegion[0].pShowText->setPosition(vPosition);	
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[1].stcOutputRegion[1].pShowText->setPosition(vPosition);
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[1].stcOutputRegion[2].pShowText->setPosition(vPosition);
		vPosition += osg::Vec3(0,-100,0);
		clsHudInstance[1].stcOutputRegion[3].pShowText->setPosition(vPosition);
	}
	pHudGroup->setUpdateCallback(new HUDUpdate);
	pRoot->addChild(pHudGroup.get());
}