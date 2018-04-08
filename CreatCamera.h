#ifndef CREATCAMERA_H
#define CREATCAMERA_H
#include "Satellite.h"
#include <osgViewer/Viewer>

void RegisterSatelliteCameraUpdate(MainBody *pSatellite);
void CreateCamera(osgViewer::Viewer& viewer);

#endif