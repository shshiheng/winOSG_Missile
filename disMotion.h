#ifndef __DISTRIBUTE_MOTION_H
#define __DISTRIBUTE_MOTION_H

#include "UserMechanics.h"

class CMotionAgent
{
public:
	SSimulationDateTime	stcGMT;
	long long			llProcessTime;
	long long			llStopTime;
	unsigned long		ulTimeStep;
	unsigned long		ulSynchronousTime;
	unsigned short		usOutputTotalNo;
	bool				bOut;
private:
	unsigned short		usCurrentNo;
public:
	void	Step(void);
	inline double		GetTimeNow(void){	return (this->llProcessTime /1000.0); };
	inline double		GetTimeStep(void){	return (this->ulTimeStep/1000000.0); };
	bool	Reset( const signed char bIsAgent, const char chrFileName[] );
	unsigned long ShiftTimeStep(signed char scRatio);
};

#endif
