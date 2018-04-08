#include <osg/LineWidth>
#include "CreateDevice.h"
#include "EarthConstant.h"

static unsigned short usFormationOrbitIndex = 0;
#define MAX_DATATYPE_NUMBER 50
//extern osg::ref_ptr<osg::Vec3Array> pOrbitVertices;//顶点
//extern osg::ref_ptr<osg::Group> pViewerRoot;
osg::Node* pCreateEarth()
{
	osg::ref_ptr<osg::MatrixTransform> matrixEarth = new osg::MatrixTransform;
	osg::ref_ptr<osg::TessellationHints> pHints = new osg::TessellationHints;
	pHints->setDetailRatio(5.0f);
	osg::ref_ptr<osg::Sphere> sphereEarth = new osg::Sphere(osg::Vec3(0.0,0.0,0.0),EarthRadius);
	osg::ref_ptr<osg::ShapeDrawable> pShapeOfEarth = new osg::ShapeDrawable(sphereEarth,pHints);
	pShapeOfEarth->setName("Earth");
	osg::Geode* pGeodeOfEarth = new osg::Geode;
	pGeodeOfEarth->addDrawable(pShapeOfEarth.get());
	//贴地图
	std::string strFileName = osgDB::findDataFile("Earth.jpg");
	osg::ref_ptr<osg::Texture2D> textureEarth = new osg::Texture2D(osgDB::readImageFile(strFileName));
	pGeodeOfEarth->getOrCreateStateSet()->setTextureAttributeAndModes(0,textureEarth);
   /* pGeodeOfEarth->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);*/
	return pGeodeOfEarth;
}

void pCreateEquator(const char chrMessages[],const unsigned char ucLength)
{
	osg::ref_ptr<osg::Geode> pGeodeOfEquator = new osg::Geode();
	osg::ref_ptr<osg::Geometry> pSurfaceOfEquator = new osg::Geometry();
	osg::ref_ptr<osg::Geometry> pCurveOfEquator = new osg::Geometry();
	osg::ref_ptr<osg::Vec3Array> pOrbitVertices = new osg::Vec3Array(36000);//顶点
	{
		pSurfaceOfEquator->setVertexArray(pOrbitVertices.get());
		osg::ref_ptr<osg::Vec4Array> pColors = new osg::Vec4Array;
		pColors->push_back(osg::Vec4(1.0f,1.0f,1.0f,0.6f));
		pSurfaceOfEquator->setColorArray(pColors.get());
		pSurfaceOfEquator->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色
		osg::ref_ptr<osg::Vec3Array> pNormals = new osg::Vec3Array;
		pNormals->push_back(osg::Vec3(0.0,-1.0,0.0));//0.0,-1.0f,0.0f
		pSurfaceOfEquator->setNormalArray(pNormals.get());
		pSurfaceOfEquator->setNormalBinding(osg::Geometry::BIND_OVERALL);//绑定法线
		pSurfaceOfEquator->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POLYGON,0,360));

		osg::StateSet* pStateSet = new osg::StateSet();
		osg::PolygonStipple* pPolygonStipple = new osg::PolygonStipple;
		pStateSet->setAttributeAndModes(pPolygonStipple,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
		osg::BlendEquation* pBlendEquation = new osg::BlendEquation(osg::BlendEquation::FUNC_ADD);    
		pStateSet->setAttributeAndModes(pBlendEquation,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
		pStateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);         
		pStateSet->setMode( GL_BLEND, osg::StateAttribute::ON );
		pSurfaceOfEquator->setStateSet(pStateSet);
	}
	{
		pCurveOfEquator->setVertexArray(pOrbitVertices.get());
		osg::Vec4Array* pColors = new osg::Vec4Array;
		pColors->push_back(osg::Vec4(1.0f,1.0f,1.0f,1.0f));
		pCurveOfEquator->setColorArray(pColors);
		pCurveOfEquator->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色
		osg::Vec3Array* pNormals = new osg::Vec3Array;
		pNormals->push_back(osg::Vec3(0.0,-1.0f,0.0f));
		pCurveOfEquator->setNormalArray(pNormals);
		pCurveOfEquator->setNormalBinding(osg::Geometry::BIND_OVERALL);//绑定法线
		pCurveOfEquator->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,360));	
	}
	pGeodeOfEquator->addDrawable(pCurveOfEquator.get());
	pGeodeOfEquator->addDrawable(pSurfaceOfEquator.get());
	pGeodeOfEquator->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF );

}
//osg::Node* pCreateEquator()
//{
//	osg::Geode* pGeodeOfequator = new osg::Geode();
//	osg::ref_ptr<osg::Geometry> pSurfaceOfEquator = new osg::Geometry();
//	osg::ref_ptr<osg::Geometry> pCurveOfEquator = new osg::Geometry();
//	osg::ref_ptr<osg::Vec3Array> pOrbitVertices = new osg::Vec3Array(36000);//顶点

//		(*pOrbitVertices)[i].set(SatellitePosition_X,SatellitePosition_Y,SatellitePosition_X);

//	{
//		pSurfaceOfEquator->setVertexArray(pOrbitVertices.get());
//		osg::ref_ptr<osg::Vec4Array> pColors = new osg::Vec4Array;
//		pColors->push_back(osg::Vec4(1.0f,1.0f,1.0f,0.6f));
//		pSurfaceOfEquator->setColorArray(pColors.get());
//		pSurfaceOfEquator->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色
//		osg::ref_ptr<osg::Vec3Array> pNormals = new osg::Vec3Array;
//		pNormals->push_back(osg::Vec3(0.0,-1.0,0.0));//0.0,-1.0f,0.0f
//		pSurfaceOfEquator->setNormalArray(pNormals.get());
//		pSurfaceOfEquator->setNormalBinding(osg::Geometry::BIND_OVERALL);//绑定法线
//		pSurfaceOfEquator->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POLYGON,0,360));
//
//		osg::StateSet* pStateSet = new osg::StateSet();
//		osg::PolygonStipple* pPolygonStipple = new osg::PolygonStipple;
//		pStateSet->setAttributeAndModes(pPolygonStipple,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
//		osg::BlendEquation* pBlendEquation = new osg::BlendEquation(osg::BlendEquation::FUNC_ADD);    
//		pStateSet->setAttributeAndModes(pBlendEquation,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
//		pStateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);         
//		pStateSet->setMode( GL_BLEND, osg::StateAttribute::ON );
//		pSurfaceOfEquator->setStateSet(pStateSet);
//	}
//	{
//		pCurveOfEquator->setVertexArray(pOrbitVertices.get());
//		osg::Vec4Array* pColors = new osg::Vec4Array;
//		pColors->push_back(osg::Vec4(1.0f,1.0f,1.0f,1.0f));
//		pCurveOfEquator->setColorArray(pColors);
//		pCurveOfEquator->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色
//		osg::Vec3Array* pNormals = new osg::Vec3Array;
//		pNormals->push_back(osg::Vec3(0.0,-1.0f,0.0f));
//		pCurveOfEquator->setNormalArray(pNormals);
//		pCurveOfEquator->setNormalBinding(osg::Geometry::BIND_OVERALL);//绑定法线
//		pCurveOfEquator->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,360));	
//	}
//	pGeodeOfequator->addDrawable(pCurveOfEquator.get());
//	pGeodeOfequator->addDrawable(pSurfaceOfEquator.get());
//	pGeodeOfequator->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
//	return pGeodeOfequator;
//}

osg::Node* CreateBaseline(osg::Vec3 StartPoint,osg::Vec3 EndPoint,osg::Geometry *pGeometry)
{
	osg::Geode* pGeodeBaseline = new osg::Geode();
	osg::ref_ptr<osg::Vec3Array> pVerticles = new osg::Vec3Array;

	pVerticles->push_back(StartPoint);
	pVerticles->push_back(EndPoint);

	pGeometry->setVertexArray(pVerticles);
	osg::ref_ptr<osg::Vec4Array> pColor = new osg::Vec4Array;
	pColor->push_back(osg::Vec4(1.0,1.0,1.0,1.0));
	pGeometry->setColorArray(pColor);
	pGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);

	/******设置线宽******/
	osg::ref_ptr<osg::LineWidth> pLineSize = new osg::LineWidth;
	pLineSize->setWidth(3);
	pGeometry->getOrCreateStateSet()->setAttributeAndModes(pLineSize.get(),osg::StateAttribute::ON);

	osg::ref_ptr<osg::Vec3Array> pNormals = new osg::Vec3Array;
	pNormals->push_back(osg::Vec3(0.0,-1.0f,0.0f));
	pGeometry->setNormalArray(pNormals);
	pGeometry->setNormalBinding(osg::Geometry::BIND_OVERALL);
	osg::ref_ptr<osg::DrawArrays> drawGeometry = new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP,0,2);
	pGeometry->addPrimitiveSet(drawGeometry);

	pGeodeBaseline->addDrawable(pGeometry);
	return pGeodeBaseline;
}

osg::Group* pCreateLight()
{
	osg::Group* pRoot = new osg::Group;	
	osg::ref_ptr<osg::StateSet> rootStateSet = new osg::StateSet;
	rootStateSet->setMode( GL_LIGHTING, osg::StateAttribute::ON );
	rootStateSet->setMode( GL_LIGHT0, osg::StateAttribute::ON );
	rootStateSet->setMode( GL_NORMALIZE, osg::StateAttribute::ON );
	osg::ref_ptr<osg::Light> myLight0 = new osg::Light;
	myLight0->setLightNum(0);
	myLight0->setPosition(osg::Vec4(1.0f,0.0f,0.0,0.0));
	myLight0->setAmbient(osg::Vec4(1.0,1.0,1.0,1.0f));
    myLight0->setDiffuse(osg::Vec4(1.0,1.0,1.0,1.0f));
	myLight0->setSpecular( osg::Vec4(1.0, 1.0, 1.0, 1) );
	osg::ref_ptr<osg::LightSource> lightS0 = new osg::LightSource;       
	lightS0->setLight(myLight0);
	pRoot->addChild(lightS0);
	return pRoot;
}

osg::StateSet* pCreate1DTextureStateToDecorate(osg::Node* pLoadedModel)
{
	const osg::BoundingSphere& BoundingSphere = pLoadedModel->getBound();

	osg::Image* pImage = new osg::Image;

	int iNumOfPixels = 1024;

	// allocate the image data, noPixels x 1 x 1 with 4 rgba floats - equivalent to a Vec4!
	pImage->allocateImage(iNumOfPixels,1,1,GL_RGBA,GL_FLOAT);
	pImage->setInternalTextureFormat(GL_RGBA);

	typedef std::vector<osg::Vec4> ColorBands;
	ColorBands MyColorbands;
	if (pLoadedModel->getName() == "SpaceCircleConfiguration")
	{
		MyColorbands.push_back(osg::Vec4(1.0f,0.0f,0.0f,0.6f));
	}
	if (pLoadedModel->getName() == "CircleUnderSatellite")
	{
		MyColorbands.push_back(osg::Vec4(0.0f,1.0f,0.0f,0.6f));
	}
	if (pLoadedModel->getName() != "SpaceCircleConfiguration" && pLoadedModel->getName() != "CircleUnderSatellite" && pLoadedModel->getName() != "MainOrbit")
	{
		MyColorbands.push_back(osg::Vec4(1.0f,0.0f,0.0f,1.0f));//
	}

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

	pStateSet->setTextureAttribute(0,pTexture,osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_1D,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_2D,osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_3D,osg::StateAttribute::OFF|osg::StateAttribute::OVERRIDE);

	pStateSet->setTextureAttribute(0,pTexgen,osg::StateAttribute::OVERRIDE);
	pStateSet->setTextureMode(0,GL_TEXTURE_GEN_S,osg::StateAttribute::ON|osg::StateAttribute::OVERRIDE);

	pStateSet->setAttribute(pMaterial,osg::StateAttribute::OVERRIDE);

	return pStateSet;
}

void OrbitUpdate(osg::Vec3 InVertices,osg::Drawable* drawable)
{
	static unsigned long ulCount = 2;
	osg::Geometry* geom = dynamic_cast<osg::Geometry*>( drawable );
	if ( !geom ) 
		return;
	osg::Vec3Array* vertices = dynamic_cast<osg::Vec3Array*>( geom->getVertexArray() );
	static unsigned long VerticesNo;
	if ( vertices )
	{
		osg::Vec3Array::iterator itr = vertices->begin();
		itr = VerticesNo+itr;
		//if(VerticesNo < 300000)
		(*itr).set(InVertices);	
		VerticesNo += 1;
		//if(VerticesNo>2700) 
		//	VerticesNo = 0;	
		drawable->dirtyBound();
		drawable->computeBound();
		vertices->dirty();
	}
}

void RelativeOrbitUpdate(osg::Vec3 InVertices,osg::Drawable* drawable)
{
	static unsigned long ulCount = 2;
	osg::Geometry* geom = dynamic_cast<osg::Geometry*>( drawable );
	if ( !geom ) 
		return;
	osg::Vec3Array* vertices = dynamic_cast<osg::Vec3Array*>( geom->getVertexArray() );
	static unsigned long VerticesNo;
	if ( vertices )
	{
		osg::Vec3Array::iterator itr = vertices->begin();
		itr = VerticesNo+itr;
		if(VerticesNo < 5000)
			(*itr).set(InVertices);	
		VerticesNo += 1;
		if(VerticesNo>4500) 
			VerticesNo = 0;	
		drawable->dirtyBound();
		drawable->computeBound();
		vertices->dirty();
	}
}
osg::Node* pCreateMainOrbit()
{
	osg::Geode* pGeodeOfOrbit = new osg::Geode();
	osg::ref_ptr<osg::Geometry> pGeometryOfOrbit = new osg::Geometry();
	osg::ref_ptr<osg::Geometry> geometryOrbitPlane = new osg::Geometry();

	osg::ref_ptr<osg::Vec3Array> pOrbitVertices = new osg::Vec3Array(10000000);//顶点
	for ( unsigned int i=0; i<10000000; ++i )
		(*pOrbitVertices)[i].set(0,0,0);

	
	geometryOrbitPlane->setVertexArray(pOrbitVertices.get());
	pGeometryOfOrbit->setVertexArray(pOrbitVertices.get());

	osg::ref_ptr<osg::Vec4Array> pColors = new osg::Vec4Array;
	pColors->push_back(osg::Vec4(1.0f,1.0f,1.0f,0.2f));//红色//1.0f,1.0f,0.0f,0.0f//黑色
	geometryOrbitPlane->setColorArray(pColors.get());
	geometryOrbitPlane->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色
	osg::ref_ptr<osg::Vec4Array> colorsOrbit = new osg::Vec4Array;
	colorsOrbit->push_back(osg::Vec4(1.0f,0.0f,0.0f,1.0f));
	pGeometryOfOrbit->setColorArray(pColors.get());
	pGeometryOfOrbit->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色

	osg::ref_ptr<osg::Vec3Array> pNormals = new osg::Vec3Array;
	pNormals->push_back(osg::Vec3(-1.0,0.0,0.0));//0.0,-1.0f,0.0f
	geometryOrbitPlane->setNormalArray(pNormals.get());
	geometryOrbitPlane->setNormalBinding(osg::Geometry::BIND_OVERALL);//绑定法线
    
	
	
 
    geometryOrbitPlane->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS,0,6000));
	pGeometryOfOrbit->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS,0,6000));

	osg::ref_ptr<osg::StateSet> pStateset = new osg::StateSet;
	osg::ref_ptr<osg::BlendEquation> pBlendEquation = new osg::BlendEquation(osg::BlendEquation::FUNC_ADD);    
	pStateset->setAttributeAndModes(pBlendEquation,osg::StateAttribute::OVERRIDE|osg::StateAttribute::ON);
	pStateset->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);         
	pStateset->setMode( GL_BLEND, osg::StateAttribute::ON );

	geometryOrbitPlane->setStateSet(pStateset.get());


	pGeometryOfOrbit->setUseDisplayList( false );
    pGeometryOfOrbit->setUseVertexBufferObjects( true );
	geometryOrbitPlane->setUseDisplayList( false );
    geometryOrbitPlane->setUseVertexBufferObjects( true );
	pGeometryOfOrbit->setInitialBound( osg::BoundingBox(osg::Vec3(-100000000.0,-100000000.0,-100000000.0), osg::Vec3(100000000.0,100000000.0,100000000.0)) );  //InitialBound
	geometryOrbitPlane->setInitialBound(osg::BoundingBox(osg::Vec3(-100000000.0,-100000000.0,-100000000.0), osg::Vec3(100000000.0,100000000.0,100000000.0)));
	pGeodeOfOrbit->addDrawable(pGeometryOfOrbit.get());
	pGeodeOfOrbit->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
	return pGeodeOfOrbit;
}

osg::Node* pCreateDynamicOrbit()
{
	osg::Geode* pGeodeOfOrbit = new osg::Geode();
	osg::ref_ptr<osg::Geometry> pGeometryOfOrbit = new osg::Geometry();

	osg::ref_ptr<osg::Vec3Array> pOrbitVertices = new osg::Vec3Array(2);
	(*pOrbitVertices)[0].set(0,0,0);
	(*pOrbitVertices)[1].set(1,0,0);

	osg::ref_ptr<osg::Vec4Array> pColors = new osg::Vec4Array;
	pColors->push_back(osg::Vec4(1.0f,1.0f,1.0f,0.2f));//红色//1.0f,1.0f,0.0f,0.0f//黑色
	pGeometryOfOrbit->setVertexArray(pOrbitVertices.get());
	osg::ref_ptr<osg::Vec4Array> colorsOrbit = new osg::Vec4Array;
	colorsOrbit->push_back(osg::Vec4(1.0f,0.0f,0.0f,1.0f));
	pGeometryOfOrbit->setColorArray(pColors.get());
	pGeometryOfOrbit->setColorBinding(osg::Geometry::BIND_OVERALL);//绑定颜色

	pGeometryOfOrbit->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2));

	pGeometryOfOrbit->setUseDisplayList( false );
    pGeometryOfOrbit->setUseVertexBufferObjects( true );
	pGeometryOfOrbit->setInitialBound( osg::BoundingBox(osg::Vec3(-100000000.0,-100000000.0,-100000000.0), osg::Vec3(100000000.0,100000000.0,100000000.0)) );  //InitialBound
	pGeodeOfOrbit->addDrawable(pGeometryOfOrbit.get());
	pGeodeOfOrbit->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF );

	return pGeodeOfOrbit;
}

void CreateProjectionMatrix(osg::Camera *camera)
{
	double left,right,bottom,top,zNear,zFar;
	left=-640;
	right=640;
	bottom=-512;
	top=512;
	zNear=1000000;
	zFar=100000000;
	camera->setProjectionMatrixAsFrustum(left*100,right*100,bottom*100,top*100,8*zNear,zFar);
}

