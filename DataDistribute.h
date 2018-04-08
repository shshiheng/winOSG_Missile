//File DataDistribute contains the definitions of DataDistribute functions
//File created @ 2009-10-21

#ifndef __TEMPSATELLITE_H
#define __TEMPSATELLITE_H

#include "Satellite.h"
#include "Manipulator.h"
#include "Hand.h"
#include <math.h>

void RegisterManipulatorUpdate(SceneManipulator *currentManipulator);
void RegisterHandUpdate(Hand *pHand);
void RegisterSatelliteUpdate(MainBody *pCurrentSatellite);

unsigned char SatelliteDataUpdate(const char chrMessages[],const unsigned char ucLength);

void ManipulatorDataUpdate(const char chrMessages[],const unsigned char ucLength);
void HandDataUpdate(const char chrMessages[],const unsigned char ucLength);
void PickResponse(unsigned char ucDeviceType,unsigned char ucMainBodyIndex,unsigned char ucDeviceIndex);
void CollectSatelliteData(void);
void CollectSatelliteRTIData(void);
double CalculateWingAttitude(const double  x[3]);

#endif
