//File Satellite contains the realization of Satellite functions
//File created @ 2009-11-6

#include <stdio.h>

#include <osgDB/ReadFile>
#include <osgDB/ReadFile>
#include <osg/State>
#include <osg/Material>
#include <osg/Texture1D>
#include <osg/TexGen>
#include <osg/BlendFunc>
#include <osg/StateSet>
#include <osg/BlendEquation>
#include <osg/PolygonOffset>
#include "Satellite.h"
#include "Manipulator.h"
#include "DataDistribute.h"
#include "HLBToFixXYZ.h"
#include "disEarth.h"
#include "CreateDevice.h"
#include "UnitConstant.h"
#define MAX_MAINBODY_NUMBER 5

osgViewer::Viewer *pCurrentViewer = NULL;
static unsigned char ucMainBodyCount = 0;
MainBody *pMainBodyInstance[MAX_MAINBODY_NUMBER] = {NULL}; 
osg::ref_ptr<osg::StateSet> pInitialState = new osg::StateSet;
osg::ref_ptr<osg::StateSet> pTexture = new osg::StateSet;


double RadarEarthFixedPositionXYZ[2][3] = {0};
struct SGeographyElements EarthRadarGeographyElements[2] = {{0,0,0},{0,0,0}};
double MaximumDetectDistance[2] = {0};
double MinimumElevationAngle[2] = {0};

Device::Device()
{
	this->bVisibleOfTemperature = true;
	this->bHighLight = false;
	this->bVisible = true;
	this->pMatrixTransform = new osg::MatrixTransform;
	this->pTransformMatrix = new osg::MatrixTransform;
	this->pNodeModel = new osg::Node;
	this->pNodeToShowTemperature = new osg::Geode;
	this->pTextOfTemperature = new osgText::Text;
	this->strDiscription = "";
	this->TCHRNodelsName[0] = L'\0';
	this->TCHRDeviceUnitage[0] = L'\0';
	this->DeviceParameter = 0;
	

	this->pNodeToShowTemperature->addDrawable(this->pTextOfTemperature);
	this->pMatrixTransform->addChild(this->pNodeToShowTemperature);
	pTextOfTemperature->setFont("fonts/simhei.ttf");
	pTextOfTemperature->setCharacterSize(20.0f);
	pTextOfTemperature->setAlignment(osgText::Text::CENTER_TOP);
	pTextOfTemperature->setAxisAlignment(osgText::Text::SCREEN);
	pTextOfTemperature->setColor(osg::Vec4(1.0,0.0,0.5,1.0));
	//pTextOfTemperature->setCharacterSizeMode(osgText::Text::SCREEN_COORDS);

}
osg::StateSet* pCreate1DTextureStateToDecorate(osg::Node* pLoadedModel,osg::Vec4 vTemperatureColor)
{
	const osg::BoundingSphere& BoundingSphere = 
		pLoadedModel->getBound();

	osg::Image* pImage = new osg::Image;

	int iNumOfPixels = 1024;

	// allocate the image data, noPixels x 1 x 1 with 4 rgba floats - equivalent to a Vec4!
	pImage->allocateImage(iNumOfPixels,1,1,GL_RGBA,GL_FLOAT);
	pImage->setInternalTextureFormat(GL_RGBA);

	typedef std::vector<osg::Vec4> ColorBands;
	ColorBands MyColorbands;
	MyColorbands.push_back(vTemperatureColor);


    

	double Nobands = MyColorbands.size();
	double Delta = Nobands/(float)iNumOfPixels;
	double Position = 0.0f;

	// fill in the image data.    
	osg::Vec4* pDataPtr = (osg::Vec4*)pImage->data();
	for(int iCount = 0;iCount < iNumOfPixels ; ++iCount, Position += Delta)
	{
		osg::Vec4 MyColor = MyColorbands[(int)Position];
		*pDataPtr++ = MyColor;
	}

	osg::Texture1D* pTexture = new osg::Texture1D;
	pTexture->setWrap(osg::Texture1D::WRAP_S,osg::Texture1D::MIRROR);
	pTexture->setFilter(osg::Texture1D::MIN_FILTER,osg::Texture1D::LINEAR);
	pTexture->setImage(pImage);

	double ZBase = BoundingSphere.center().z()-BoundingSphere.radius();
	double ZScale = 2.0f/BoundingSphere.radius();

	osg::TexGen* pTexgen = new osg::TexGen;
	pTexgen->setMode(osg::TexGen::OBJECT_LINEAR);
	pTexgen->setPlane(osg::TexGen::S,osg::Plane(0.0f,0.0f,ZScale,-ZBase));

	osg::Material* pMaterial = new osg::Material;

	osg::StateSet* pStateSet = new osg::StateSet;
	pStateSet->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
	pStateSet->setTextureAttribute(0,pTexture,osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_1D,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_2D,osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_3D,osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);

	pStateSet->setTextureAttribute(0,pTexgen,osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_GEN_S,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);

	pStateSet->setAttribute(pMaterial,osg::StateAttribute::OVERRIDE);
	//设置其透明度
	osg::BlendEquation* pBlendEquation = new osg::BlendEquation(osg::BlendEquation::FUNC_ADD);    

    pStateSet->setAttributeAndModes(pBlendEquation,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
    pStateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);         
    pStateSet->setMode( GL_BLEND, osg::StateAttribute::ON );

	return pStateSet;
}
signed char CheckFireStatus()
{
	unsigned char ucCount = 0;
	std::string strCompare;
	osg::ref_ptr<osg::StateSet> pOnFireState = new osg::StateSet;
	osg::ref_ptr<osg::StateSet> pCheckState = new osg::StateSet;
	static osg::ref_ptr<osg::StateSet> pInitialState = new osg::StateSet;

	for(ucCount = 0;ucCount<pMainBodyInstance[0]->ucActiveDeviceNumber;ucCount++)
	{
		pCheckState = pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform->getOrCreateStateSet();
		if (pCheckState.get() == pInitialState.get()&&pMainBodyInstance[0]->clsActiveDevices[ucCount].bThrusterOnFire == false)	continue;

		pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform->setDataVariance(Object::DataVariance::DYNAMIC);
		pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform->getOrCreateStateSet()->setDataVariance(Object::DataVariance::DYNAMIC);
		strCompare=pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform->getName().substr(0,8);

		if (strCompare == "Thruster"&&pMainBodyInstance[0]->clsActiveDevices[ucCount].bThrusterOnFire == true)
		{
			pOnFireState = pCreate1DTextureStateToDecorate(pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform,osg::Vec4d(255,0,0,0.5));
			pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform->setStateSet(pOnFireState.get());			
		}
		if (strCompare == "Thruster"&&pMainBodyInstance[0]->clsActiveDevices[ucCount].bThrusterOnFire == false)
			pMainBodyInstance[0]->clsActiveDevices[ucCount].pMatrixTransform->setStateSet(pInitialState.get());
	}

	return 0;
}
void Device::HighLight()
{
	pTexture = pCreate1DTextureStateToDecorate(this->pNodeModel,osg::Vec4(1.0,1.0,1.0,0.3));
	this->pMatrixTransform->setStateSet(pTexture.get());
	this->bHighLight = true;
}

void Device::CancelHighLight()
{
	this->pMatrixTransform->setStateSet(pInitialState.get());
	this->bHighLight = false;
}


void Device::DisplayTemperature(const double Temperature)
{

	osg::Vec4 vColor;
	if(Temperature<=-100)		
		vColor.set(0.0,0.0,1.0,0.3);//  <=-100
	if(Temperature>-100 &&Temperature<=-50)
		vColor.set(0.0,0.0028*Temperature+0.28,1.0,0.3);//  <=-100   0.14
	if(Temperature>-50 &&Temperature<=-25)
		vColor.set(0.0,0.0056*Temperature+0.42,1.0,0.3);  //    0
	if(Temperature>-25 &&Temperature<=0)
		vColor.set(0.0,0.0056*Temperature+0.42,1.0,0.3);  //    0
	if(Temperature>0 &&Temperature<=5)
		vColor.set(0.0,0.028*Temperature+0.42,1.0,0.3);  //    0
	if(Temperature>5 &&Temperature<=10)
		vColor.set(0.0,0.028*Temperature+0.42,1.0,0.3);  //    0
	if(Temperature>10 &&Temperature<=15)
		vColor.set(0.0,0.028*Temperature+0.42,1.0,0.3);  //    0
	if(Temperature>15 &&Temperature<=20)
		vColor.set(0.0,0.032*Temperature+0.36,1.0,0.3);  //    0
	if(Temperature>20 &&Temperature<=25)
		vColor.set(0.0,1.0,-0.068*Temperature+2.36,0.3);  //    0
	if(Temperature>25 &&Temperature<=30)
		vColor.set(0.0,1.0,-0.066*Temperature+2.31,0.3);  //    0
	if(Temperature>30 &&Temperature<=35)
		vColor.set(0.0,1.0,-0.066*Temperature+2.31,0.3);  //    0

	if(Temperature>35 &&Temperature<=40)
		vColor.set(0.05*Temperature-1.75,1.0,0.0,0.3);  //    0
	if(Temperature>40 &&Temperature<=45)
		vColor.set(0.05*Temperature-1.75,1.0,0.0,0.3);  //    0
	if(Temperature>45 &&Temperature<=50)
		vColor.set(0.05*Temperature-1.75,1.0,0.0,0.3);  //    0
	if(Temperature>50 &&Temperature<=55)
		vColor.set(0.05*Temperature-1.75,1.0,0.0,0.3);  //    0
	if(Temperature>55 &&Temperature<=60)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>60 &&Temperature<=65)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>65 &&Temperature<=70)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>70 &&Temperature<=75)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>75 &&Temperature<=80)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>80 &&Temperature<=90)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>90 &&Temperature<=100)
		vColor.set(1.0,-0.025*Temperature+2.375,0.0,0.3);  //    0
	if(Temperature>100)
		vColor.set(1.0,0.0,0.0,0.3);  //    0
	pTexture = pCreate1DTextureStateToDecorate(this->pNodeModel,vColor);
	this->pNodeModel->setStateSet(pTexture.get());
}

DeviceWithAction::DeviceWithAction()
{
	this->vPosition = osg::Vec3(0,0,0);
	this->vEulerAttitude = osg::Vec3(0,0,0);
	this->vScale = osg::Vec3(1,1,1);
	this->mFixMatrix = osg::Matrix::identity();
	this->qQuaternionAttitude = osg::Quat(0,0,0,1);
	this->vSpeed=osg::Vec3(0,0,0);
	this->vAngularSpeed = osg::Vec3(0,0,0);
	this->bRotateByQuaternion = true;
	this->CurrentDeviceParameter = 0;
	this->TCHRCurrentDeviceName[0]  =L'\0'; 
	this->bThrusterOnFire = false;
}
void DeviceWithAction::operator ()(osg::Node* node,osg::NodeVisitor* nv)
{
	this->vEulerAttitude += this->vAngularSpeed * FRAME_TIME_RADIO;
	this->pTransformMatrix = dynamic_cast<osg::MatrixTransform*>(node);
	osg::Matrix mTranslate,mRotate,mScale;
	mTranslate.makeTranslate(this->vPosition);
	if(this->bRotateByQuaternion)
		mRotate.makeRotate(this->qQuaternionAttitude);
	else
	mRotate.makeRotate(osg::DegreesToRadians(this->vEulerAttitude._v[2]),osg::Vec3(0,0,1),osg::DegreesToRadians(this->vEulerAttitude._v[0]),osg::Vec3(1,0,0),osg::DegreesToRadians(this->vEulerAttitude._v[1]),osg::Vec3(0,1,0));

	mScale.makeScale(this->vScale);
	pTransformMatrix->setMatrix(mRotate*mScale*(this->mFixMatrix)*mTranslate);
	traverse(node,nv);
}

MainBody::MainBody( osg::Group *pRoot,osgViewer::Viewer *pViewer)
{
	pCurrentViewer = pViewer;
	this->ucDeviceNumber = 0;
	this->ucActiveDeviceNumber = 0;
	pMainBodyInstance[ucMainBodyCount] = this;
	ucMainBodyCount ++;
}

bool MainBody::ReadConfigFile(const char chrFileName[])
{
	FILE *pFileStream = fopen(chrFileName,"r");
	if(pFileStream == NULL)
		return false;
	char chrDeviceType[16]    = {0};
	char chrModelName[32] = {0};
	int		iCount		  = 0;
	int		iFatherNode	  = 0;
	int		iDisplay	  = 0;
	char chrNodelName[32] = {0};
	double FixPositionX   = 0;
	double FixPositionY   = 0;
	double FixPositionZ   = 0;
	double FixQuat0		  = 0.0;
	double FixQuat1		  = 0.0;
	double FixQuat2		  = 0.0;
	double FixQuat3		  = 0.0;
	double Scale		  = 0;
	char chrDiscription[32] = {0};
	char chrFilePath[32] = {0};
	RadarGeographyElements EarthRadar[2];
	CEarthAgent clsEarthSatellite;
	EarthRadar[0].readEarthRadarFile("EarthRadar1.INI");
	EarthRadar[1].readEarthRadarFile("EarthRadar3.INI");
	EarthRadarGeographyElements[0].Longitude = EarthRadar[0].Longitude/180*osg::PI;
	EarthRadarGeographyElements[0].Latitude = EarthRadar[0].Latitude/180*osg::PI;
	EarthRadarGeographyElements[0].Altitude = EarthRadar[0].Altitude*ToM;
	EarthRadarGeographyElements[1].Longitude = EarthRadar[1].Longitude/180*osg::PI;
	EarthRadarGeographyElements[1].Latitude = EarthRadar[1].Latitude/180*osg::PI;
	EarthRadarGeographyElements[1].Altitude = EarthRadar[1].Altitude*ToM;
	MaximumDetectDistance[0] = EarthRadar[0].MaximumDistance;
	MinimumElevationAngle[0] = EarthRadar[0].MinimumAngle;
	MaximumDetectDistance[1] = EarthRadar[1].MaximumDistance;
	MinimumElevationAngle[1] = EarthRadar[1].MinimumAngle;
	clsEarthSatellite.HLBToEarthFixedPosition(EarthRadarGeographyElements[0],RadarEarthFixedPositionXYZ[0]);
	clsEarthSatellite.HLBToEarthFixedPosition(EarthRadarGeographyElements[1],RadarEarthFixedPositionXYZ[1]);

	
	printf("Loading Models, Please Wait...\n");
	
	fscanf(pFileStream,"DeviceType	ModelName	Index	FatherNode	Display	 NodeName	 Position	        Attitude	Scale	 Discription	FIELD");
	while(true)
	{
		fscanf(pFileStream,"%s %s %d %d %d %s (%lf,%lf,%lf) (%lf,%lf,%lf,%lf) %lf %s",chrDeviceType,chrModelName,&iCount,&iFatherNode,&iDisplay,chrNodelName,&FixPositionX,&FixPositionY,&FixPositionZ,&FixQuat0,&FixQuat1,&FixQuat2,&FixQuat3,&Scale,chrDiscription);

		if(strcmp(chrDeviceType,"Device") == 0)
		{
			memset(chrFilePath,0,32);
			memcpy(chrFilePath, "Model/",8);
			this->clsDevices[iCount].pNodeModel = osgDB::readNodeFile(strcat(chrFilePath,chrModelName));
			this->clsDevices[iCount].pNodeModel->getOrCreateStateSet()->setMode(GL_NORMALIZE,osg::StateAttribute::ON);
			this->clsDevices[iCount].pNodeModel->setNodeMask(iDisplay);
			this->clsDevices[iCount].strDiscription += chrDiscription;
			this->clsDevices[iCount].pMatrixTransform->addChild(this->clsDevices[iCount].pNodeModel);
			this->clsDevices[iCount].pMatrixTransform->setName(chrNodelName);
			this->clsDevices[iCount].pMatrixTransform->setMatrix(osg::Matrix::rotate(osg::Quat(FixQuat1,FixQuat2,FixQuat3,FixQuat0))*osg::Matrix::scale(osg::Vec3(Scale,Scale,Scale))*osg::Matrix::translate(FixPositionX,FixPositionY,FixPositionZ));
			if(iFatherNode == -1)
				this->pMatrixTransform->addChild(this->clsDevices[iCount].pMatrixTransform);
			else 
				this->clsDevices[iFatherNode].pMatrixTransform->addChild(this->clsDevices[iCount].pMatrixTransform);
			this->ucDeviceNumber ++;
		}
		else if(strcmp(chrDeviceType,"ActiveDevice") == 0)
		{
			memset(chrFilePath,0,32);
			memcpy(chrFilePath, "Model/",8);
			this->clsActiveDevices[iCount].pNodeModel = osgDB::readNodeFile(strcat(chrFilePath,chrModelName));
			this->clsActiveDevices[iCount].pNodeModel->getOrCreateStateSet()->setMode(GL_NORMALIZE,osg::StateAttribute::ON);
			this->clsActiveDevices[iCount].pNodeModel->setNodeMask(iDisplay);
			this->clsActiveDevices[iCount].strDiscription += chrDiscription;
			this->clsActiveDevices[iCount].mFixMatrix = osg::Matrix::rotate(osg::Quat(FixQuat1,FixQuat2,FixQuat3,FixQuat0))*osg::Matrix::scale(osg::Vec3(Scale,Scale,Scale))*osg::Matrix::translate(FixPositionX,FixPositionY,FixPositionZ);
			this->clsActiveDevices[iCount].pMatrixTransform->addChild(this->clsActiveDevices[iCount].pNodeModel);
			this->clsActiveDevices[iCount].pMatrixTransform->setName(chrNodelName);
			this->clsActiveDevices[iCount].pMatrixTransform->setUpdateCallback(&(this->clsActiveDevices[iCount]));

			if(iFatherNode == -1)
				this->pMatrixTransform->addChild(this->clsActiveDevices[iCount].pMatrixTransform);
			else
				this->clsActiveDevices[iFatherNode].pMatrixTransform->addChild(this->clsActiveDevices[iCount].pMatrixTransform);
			this->ucActiveDeviceNumber ++;
			
		}
		else if(strcmp(chrDeviceType,"End") == 0)
			break;
	}
	this->clsActiveDevices[80].pNodeModel = pCreateMainOrbit();
	this->clsActiveDevices[80].pNodeModel->setNodeMask(1);
	this->clsActiveDevices[80].pMatrixTransform->addChild(this->clsActiveDevices[80].pNodeModel);
	this->pMatrixTransform->addChild(this->clsActiveDevices[80].pMatrixTransform);	

	this->clsActiveDevices[81].pNodeModel = pCreateMainOrbit();
	this->clsActiveDevices[81].pNodeModel->setNodeMask(1);
	this->clsActiveDevices[81].pMatrixTransform->addChild(this->clsActiveDevices[81].pNodeModel);
	this->pMatrixTransform->addChild(this->clsActiveDevices[81].pMatrixTransform);	

	this->clsActiveDevices[82].pNodeModel = pCreateMainOrbit();
	this->clsActiveDevices[82].pNodeModel->setNodeMask(1);
	this->clsActiveDevices[82].pMatrixTransform->addChild(this->clsActiveDevices[82].pNodeModel);
	this->clsActiveDevices[82].pMatrixTransform->setUpdateCallback(&(this->clsActiveDevices[82]));
	this->pMatrixTransform->addChild(this->clsActiveDevices[82].pMatrixTransform);

	this->clsActiveDevices[83].pNodeModel = pCreateEarth();
	this->clsActiveDevices[83].pNodeModel->setNodeMask(1);
	this->clsActiveDevices[83].pMatrixTransform->addChild(this->clsActiveDevices[83].pNodeModel);
	this->clsActiveDevices[83].pMatrixTransform->setUpdateCallback(&(this->clsActiveDevices[83]));
	this->clsActiveDevices[83].bRotateByQuaternion = false;
	this->clsActiveDevices[83].vEulerAttitude[2] = 180;
	this->pMatrixTransform->addChild(this->clsActiveDevices[83].pMatrixTransform);


	osg::ref_ptr<osg::Material> matirial = new osg::Material();
	osg::ref_ptr<osg::StateSet>EarthStateSet = new osg::StateSet();
	EarthStateSet = this->clsActiveDevices[83].pNodeModel->getOrCreateStateSet();
	matirial->setColorMode(osg::Material::DIFFUSE);
	matirial->setAmbient(osg::Material::FRONT, osg::Vec4(0.2,0.2,0.2,1));
	matirial->setDiffuse(osg::Material::FRONT, osg::Vec4(0.9,0.9,0.9,1));
	matirial->setSpecular(osg::Material::FRONT, osg::Vec4(0.5, 0.5, 0.5, 1));
	matirial->setShininess(osg::Material::FRONT, 128.0f);
	matirial->setEmission(osg::Material::FRONT,osg::Vec4(0, 0, 0, 1));
	EarthStateSet->setAttributeAndModes(matirial.get(), osg::StateAttribute::ON); 


	fclose(pFileStream);
	printf("Loading complete!\n");

	return true;
}
void MainBody::AssembleSatellite()
{
	;
}


bool Pick(double MouseX,double MouseY)
{
    osgUtil::LineSegmentIntersector::Intersections Intersections;
	std::string strOutput="";

	bool bMainBodySelected = false;
	bool bDeviceSelected = false;
	bool bActiveDeviceSelected =false;
	unsigned char ucMainBodySelected = 0;
	unsigned char ucSelectDeviceNo = 0;
	unsigned char ucSelectActiveDeviceNo = 0;

	double CenterX = MouseX;
	double CenterY = MouseY;

#ifdef DEBUG
	cout<<"PICK has been functioned"<<endl;
#endif

	if(pCurrentViewer->computeIntersections(CenterX,CenterY,Intersections))
	{
		osgUtil::LineSegmentIntersector::Intersections::iterator HitIterator = Intersections.begin();
		for(HitIterator; HitIterator != Intersections.end(); ++HitIterator )
		{
			 if (!HitIterator->nodePath.empty() && !(HitIterator->nodePath.back()->getName().empty()))
			 { 
				 cout<<HitIterator->nodePath.back()->getName()<<endl;
				 osg::NodePath nodePath = HitIterator->nodePath;

				 for(int iSize = nodePath.size() - 1;iSize >= 0;iSize--)
				 {
					 osg::Node *pNode = dynamic_cast<osg::Node *>(nodePath[iSize]);
					 cout<<pNode->getName()<<endl;
					 if(!pNode->getName().empty())
					 {
						 for(unsigned char ucIndex = 0;ucIndex <  ucMainBodyCount;ucIndex ++)
						 {
							 if(pNode->getName() == pMainBodyInstance[ucIndex]->pMatrixTransform->getName())
							 {
								 ucMainBodySelected = ucIndex;
								 bMainBodySelected = true;
								 bDeviceSelected = false;
								 bActiveDeviceSelected = false;
								 PickResponse(0,ucMainBodySelected,0);
#ifdef DEBUG
									 cout<<pNode->getName()<<"   @ "<<iSize<<" floor"<<endl;
									 cout<<"Discription: "<<pMainBodyInstance[ucIndex]->strDiscription<<endl;
#endif	
								 return true;
							 }
							 for(unsigned char ucDeviceIndex = 0;ucDeviceIndex < pMainBodyInstance[ucIndex]->ucDeviceNumber;ucDeviceIndex ++)
							 {
								 if(pNode->getName() ==  pMainBodyInstance[ucIndex]->clsDevices[ucDeviceIndex].pMatrixTransform->getName())
								 {
									 ucMainBodySelected = ucIndex;

									 bMainBodySelected = false;
									 bDeviceSelected = true;
									 bActiveDeviceSelected = false;
									 ucSelectDeviceNo = ucDeviceIndex;
									 PickResponse(1,bMainBodySelected,ucSelectDeviceNo);
									
#ifdef DEBUG
									 cout<<pNode->getName()<<"   @ "<<iSize<<" floor"<<endl;
									 cout<<"Discription: "<<pMainBodyInstance[ucIndex]->clsDevices[ucDeviceIndex].strDiscription<<endl;
#endif
									 return true;
								 }
							 }
							 for(unsigned char ucActiveDeviceIndex = 0;ucActiveDeviceIndex < pMainBodyInstance[ucIndex]->ucActiveDeviceNumber;ucActiveDeviceIndex++)
							 {
								 if(pNode->getName() == pMainBodyInstance[ucIndex]->clsActiveDevices[ucActiveDeviceIndex].pMatrixTransform->getName())
								 {
									 ucMainBodySelected = ucIndex;

									 bMainBodySelected = false;
									 bDeviceSelected = false;
									 bActiveDeviceSelected = true;
									 ucSelectActiveDeviceNo = ucActiveDeviceIndex;
									 PickResponse(2,ucMainBodySelected,ucSelectActiveDeviceNo);
#ifdef DEBUG
									 cout<<pNode->getName()<<"   @ "<<iSize<<" floor"<<endl;
									 cout<<"Discription: "<<pMainBodyInstance[ucIndex]->clsActiveDevices[ucActiveDeviceIndex].strDiscription<<endl;
#endif
									 return true;
								 }
							 }
						 }
					 }
				 }
			 }
		}
	}
	bMainBodySelected = false;
	bDeviceSelected = false;
	bActiveDeviceSelected = false;
	return false;

}

