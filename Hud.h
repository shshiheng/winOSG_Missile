//File Hud contains the definitions of Hud functions
//File created @ 2009-10-30

#ifndef __HUD_H
#define __HUD_H

#include <iostream>

#include <osg/Node>
#include <osgText/Text>
#include <osgViewer/Viewer>
using namespace std;


void AddHUDToScene(osg::Group* pRoot);
void WriteToHud(const unsigned char ucHudIndex,const unsigned char ucRegionIndex,const string strOutput);
void WriteToScrollHud(const unsigned char ucScrollHudIndex,string strShowContent);

#endif 