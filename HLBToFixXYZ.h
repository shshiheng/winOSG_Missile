#ifndef __HLBTOFIXXYZ_H
#define __HLBTOFIXXYZ_H

#include<math.h>
#include<stdio.h>
class RadarGeographyElements
{
public:
	double Longitude;
	double Latitude;
	double Altitude;
	double MaximumDistance;
	double MinimumAngle;
public:
	RadarGeographyElements();
	RadarGeographyElements(double Longitude,double Latitude,double Altitude,double MaximumDetectDistance,double MinimumElevationAngle);
	~RadarGeographyElements();
	bool readEarthRadarFile(const char chrFileName[]);
};
void StationAntennaRotateToSatellite(const  double  Frame1[3],const double  Frame2[3],double QuaternionAttitude[4],double *StationToSateAngle);

#endif