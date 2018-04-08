#include "HLBToFixXYZ.h"

void StationAntennaRotateToSatellite(const  double  Frame1[3],const double  Frame2[3],double QuaternionAttitude[4],double *StationToSateAngle)
{
	double x,y,z,r,CosB[3]; 
	x = Frame1[1]*Frame2[2]-Frame1[2]*Frame2[1];
	y = Frame1[2]*Frame2[0]-Frame1[0]*Frame2[2];
	z = Frame1[0]*Frame2[1]-Frame1[1]*Frame2[0];
	r = sqrt(x*x+y*y+z*z);
	CosB[0] = x/r;CosB[1] = y/r;CosB[2] = z/r;
	*StationToSateAngle = acos((Frame1[0]*Frame2[0]+Frame1[1]*Frame2[1]+Frame1[2]*Frame2[2])/(sqrt(Frame1[0]*Frame1[0]+Frame1[1]*Frame1[1]+Frame1[2]*Frame1[2])+sqrt(Frame2[0]*Frame2[0]+Frame2[1]*Frame2[1]+Frame2[2]*Frame2[2])));
	QuaternionAttitude[0] = cos(*StationToSateAngle/2);
	QuaternionAttitude[1] = sin(*StationToSateAngle/2)*CosB[0];
	QuaternionAttitude[2] = sin(*StationToSateAngle/2)*CosB[1];
	QuaternionAttitude[3] = sin(*StationToSateAngle/2)*CosB[2];
}

RadarGeographyElements::RadarGeographyElements()
{
	Longitude = 0;
	Latitude = 0;
	Altitude = 0;
	MaximumDistance = 0;
	MinimumAngle = 0;
}
RadarGeographyElements::RadarGeographyElements(double Longitude,double Latitude,double Altitude,double MaximumDetectDistance,double MinimumElevationAngle)
{
	Longitude = Longitude;
	Latitude = Latitude;
	Altitude = Altitude;
	MaximumDetectDistance = MaximumDetectDistance;
	MinimumElevationAngle = MinimumElevationAngle;
}
RadarGeographyElements::~RadarGeographyElements()
{
}
bool RadarGeographyElements::readEarthRadarFile(const char chrFileName[])
{
	FILE *pFileStream = fopen(chrFileName,"r");
	if(pFileStream == NULL)
		return false;
	fscanf(pFileStream,"Geography Location = %lf km %lf deg %lf deg\n",&Altitude,&Longitude,&Latitude);
	fscanf(pFileStream,"Maximum Detect Distance = %lf km\n",&MaximumDistance);
	fscanf(pFileStream,"Minimum Elevation Angle = %lf deg\n",&MinimumAngle);
	return true;
}