#include "Hand.h"

#include <stdio.h>
#include <iostream>
using namespace std;
using namespace osg;

void HandVisitor::apply(osg::Node& Node)
{
	if(Node.getName() == "thumb")		//��Ĵָ��չ״̬
	{
		this->pNodeThumb = &Node;			//�ѽ��ܽ������ݴ��͹���������֪���Ǹ�DB��㣬����ֱ�����
	}
	if(Node.getName() == "thumbFlex")//����״̬
	{
		this->pNodeThumbFlex = &Node;
	}
	if(Node.getName() == "index")		//ʳָ
	{
		this->pNodeIndex = &Node;
	}
	if(Node.getName() == "indexFlex")
	{
		this->pNodeIndexFlex = &Node;
	}
	if(Node.getName() == "middle")	//��ָ
	{
		this->pNodeMiddle = &Node;
	}
	if(Node.getName() == "middleFlex")
	{
		this->pNodeMiddleFlex = &Node;
	}
	if(Node.getName() == "ring")		//����ָ
	{
		this->pNodeRing = &Node;
	}
	if(Node.getName() == "ringFlex")
	{
		this->pNodeRingFlex = &Node;
	}
	if(Node.getName() == "little")	//Сָ
	{
		this->pNodeLittle = &Node;
	}
	if(Node.getName() == "littleFlex")
	{
		this->pNodeLittleFlex = &Node;
	}
	traverse(Node);					//���´���
}

void Hand::operator()(osg::Node* pNode,osg::NodeVisitor* pNodeVisitor)
{
	osg::MatrixTransform* m_part=dynamic_cast<osg::MatrixTransform*>(pNode);

	osg::Matrixd mTranslate,mRotate;
	mTranslate.makeTranslate(this->HandPosition);
	mRotate= osg::Matrixd::rotate(osg::DegreesToRadians(this->HandAttitude._v[2]),osg::Vec3(0,0,1),osg::DegreesToRadians(this->HandAttitude._v[0]),osg::Vec3(1,0,0),osg::DegreesToRadians(this->HandAttitude._v[1]),osg::Vec3(0,1,0));

	m_part->setMatrix(osg::Matrix::scale(30, 30, 30)*mTranslate*mRotate);

	if(this->bFingleState[0] == 0) this->clsHandNodeVisitor.pNodeThumb->setNodeMask(1);			//0��ʱ����չ״̬�Ĵ�ķָ������
	else if (this->bFingleState[0] == 1) this->clsHandNodeVisitor.pNodeThumb->setNodeMask(0);	//1��ʱ����չ״̬�Ĵ�ķָ����

	if(this->bFingleState[0] == 0) this->clsHandNodeVisitor.pNodeThumbFlex->setNodeMask(0);
	else if (this->bFingleState[0] == 1) this->clsHandNodeVisitor.pNodeThumbFlex->setNodeMask(1);

	if(this->bFingleState[1] == 0) this->clsHandNodeVisitor.pNodeIndex->setNodeMask(1);			//ʳָ
	else if (this->bFingleState[1] == 1) this->clsHandNodeVisitor.pNodeIndex->setNodeMask(0);

	if(this->bFingleState[1] == 0) this->clsHandNodeVisitor.pNodeIndexFlex->setNodeMask(0);
	else if (this->bFingleState[1] == 1) this->clsHandNodeVisitor.pNodeIndexFlex->setNodeMask(1);

	if(this->bFingleState[2] == 0) this->clsHandNodeVisitor.pNodeMiddle->setNodeMask(1);		//��ָ
	else if (this->bFingleState[2] == 1) this->clsHandNodeVisitor.pNodeMiddle->setNodeMask(0);

	if(this->bFingleState[2] == 0) this->clsHandNodeVisitor.pNodeMiddleFlex->setNodeMask(0);
	else if (this->bFingleState[2] == 1) this->clsHandNodeVisitor.pNodeMiddleFlex->setNodeMask(1);

	if(this->bFingleState[3] == 0) this->clsHandNodeVisitor.pNodeRing->setNodeMask(1);			//����ָ
	else if (this->bFingleState[3] == 1) this->clsHandNodeVisitor.pNodeRing->setNodeMask(0);

	if(this->bFingleState[3] == 0) this->clsHandNodeVisitor.pNodeRingFlex->setNodeMask(0);
	else if (this->bFingleState[3] == 1) this->clsHandNodeVisitor.pNodeRingFlex->setNodeMask(1);

	if(this->bFingleState[4] == 0) this->clsHandNodeVisitor.pNodeLittle->setNodeMask(1);		//Сָ
	else if (this->bFingleState[4] == 1) this->clsHandNodeVisitor.pNodeLittle->setNodeMask(0);

	if(this->bFingleState[4] == 0) this->clsHandNodeVisitor.pNodeLittleFlex->setNodeMask(0);
	else if (this->bFingleState[4] == 1) this->clsHandNodeVisitor.pNodeLittleFlex->setNodeMask(1);

	traverse(pNode,pNodeVisitor);
}

void cHandFollowViewpoint::operator()(osg::Node* pNode,osg::NodeVisitor* pNodeVisitor)
{
	osg::PositionAttitudeTransform* pMatrix =dynamic_cast<osg::PositionAttitudeTransform*>(pNode);

//Code here aim to seperate the hand model some distance away from view point
	double HandViewPointDistance = 100;
	osg::Vec3 vDistance;
	vDistance._v[0] = HandViewPointDistance*cosf(pCurrentManipulator->CameraLatitude)*sinf(pCurrentManipulator->CameraLongitude);
	vDistance._v[1] = -HandViewPointDistance*cosf(pCurrentManipulator->CameraLatitude)*cosf(pCurrentManipulator->CameraLongitude);
	vDistance._v[2] = HandViewPointDistance*sinf(pCurrentManipulator->CameraLatitude);

	osg::Vec3 vHandPosition = pCurrentManipulator->vViewPointPosition - vDistance;
	pMatrix->setPosition(vHandPosition);
	osg::Matrixd mCamera = pCurrentViewer->getCameraManipulator()->getMatrix();
	pCurrentViewer->getCamera()->getProjectionMatrix();

	osg::Quat qHandAttitude = mCamera.getRotate();
	//�Ӿ�����ϵ��ͶӰ����ϵY��Z�ụ������Ҫת��90��
	qHandAttitude = osg::Quat(osg::DegreesToRadians(-90.0),osg::Vec3(1,0,0)) * qHandAttitude;
	pMatrix->setAttitude(qHandAttitude);
	traverse(pNode,pNodeVisitor);
}

Hand::Hand()
{
	this->bFingleState[0] = 0;
	this->bFingleState[1] = 0;
	this->bFingleState[2] = 0;
	this->bFingleState[3] = 0;
	this->bFingleState[4] = 0;

	this->HandPosition = osg::Vec3(0,0,0);
	this->HandAttitude = osg::Vec3(0,0,0);
}
