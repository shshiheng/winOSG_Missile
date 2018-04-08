#include "CreatCamera.h"
#include "Satellite.h"

MainBody *pCurrentSatelliteCamera = NULL;
void RegisterSatelliteCameraUpdate(MainBody *pSatellite)
{
	pCurrentSatelliteCamera = pSatellite;
}
void CreateCamera(osgViewer::Viewer& viewer)
{
	osg::GraphicsContext::WindowingSystemInterface* wsi = osg::GraphicsContext::getWindowingSystemInterface();

	unsigned int width, height;
	wsi->getScreenResolution(osg::GraphicsContext::ScreenIdentifier(0), width, height);

		osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
		traits->x = 0;
		traits->y = 0;
		traits->width = width;
		traits->height = height;
		traits->windowDecoration = true;
		traits->doubleBuffer = true;
		traits->sharedContext = 0;

		osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());

		gc->setClearColor(osg::Vec4f(0.2f,0.2f,0.6f,1.0f));
		gc->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		//�����
		osg::ref_ptr<osg::Camera> cameraMaster = new osg::Camera;
		cameraMaster->setGraphicsContext(gc.get());
		//����ӿ�����
		cameraMaster->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
		// ��������VIEWER
		viewer.addSlave(cameraMaster.get());

		//����ͼ
		double aspectRatioScale = 1.0;
		osg::ref_ptr<osg::Camera> camera = new osg::Camera;
		camera->setGraphicsContext(gc.get());
		camera->setViewport(new osg::Viewport(0,0, width/4, height/4));
		/*camera->setViewMatrixAsLookAt(osg::Vec3(10000000,10000000,0),osg::Vec3(0,0,0),osg::Vec3(1,0,0));*/
		GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
		camera->setDrawBuffer(buffer);
		camera->setReadBuffer(buffer);
		viewer.addSlave(camera.get(), osg::Matrixd(), osg::Matrixd::rotate(osg::DegreesToRadians(180.0), 0.0, 1.0, 0.0));

	//��ʼ������ͼ
	osg::ref_ptr<osg::Camera> rightViewer = new osg::Camera;
	rightViewer->setGraphicsContext(gc.get());
	rightViewer->setViewport(new osg::Viewport(traits ->width/4,0, traits ->width/4, traits->height/4));
	GLenum bufferright = traits->doubleBuffer ? GL_BACK : GL_FRONT;
	rightViewer->setDrawBuffer(bufferright);
	rightViewer->setReadBuffer(bufferright);
	//rightViewer->setProjectionMatrix(osg::Matrixd::scale(-10.0f,1.0f,1.0f));
	viewer.addSlave(rightViewer.get(), osg::Matrixd(), osg::Matrixd::rotate(osg::DegreesToRadians(90.0), 0.0, 1.0, 0.0));


	//��ʼ��ǰ��ͼ
	osg::ref_ptr<osg::Camera> frontViewer = new osg::Camera;
	frontViewer->setGraphicsContext(gc.get());
	frontViewer->setViewport(new osg::Viewport((traits ->width/4)*2-1,0, traits ->width/4, traits->height/4));
	GLenum bufferfront = traits->doubleBuffer ? GL_BACK : GL_FRONT;
	frontViewer->setDrawBuffer(bufferfront);
	frontViewer->setReadBuffer(bufferfront);

	//rightViewer->setProjectionMatrix(osg::Matrixd::scale(-10.0f,1.0f,1.0f));
	viewer.addSlave(frontViewer.get(), osg::Matrixd(), osg::Matrixd::rotate(osg::DegreesToRadians(0.0), 0.0, 1.0, 0.0));


	//��ʼ������ͼ
	osg::ref_ptr<osg::Camera> leftViewer = new osg::Camera;
	leftViewer->setGraphicsContext(gc.get());
	leftViewer->setViewport(new osg::Viewport((traits ->width/4)*3-1,0.0, traits ->width/4, traits->height/4));
	GLenum bufferleft = traits->doubleBuffer ? GL_BACK : GL_FRONT;
	leftViewer->setDrawBuffer(bufferleft);
	leftViewer->setReadBuffer(bufferleft);
	//rightViewer->setProjectionMatrix(osg::Matrixd::scale(-10.0f,1.0f,1.0f));
	viewer.addSlave(leftViewer.get(), osg::Matrixd(), osg::Matrixd::rotate(osg::DegreesToRadians(270.0), 0.0, 1.0, 0.0));




}