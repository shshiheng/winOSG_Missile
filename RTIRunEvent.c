#include <string.h>
#include <windows.h>
#include "RTIRunEvent.h"
#include "RTI_NET.h"
#include "RTIBUSChannelNo.h"
#include "RTIHeadDefine.h"

static unsigned char	ucSynchronisedNo[MAX_STEP_NO] = {0};
static HANDLE			hRTIEvent;
static unsigned char	ucRTIFederalNo = 0;
static signed char		scIgnoreAttitude = 0;

signed char IsIgnoreAttitude()
{
	return scIgnoreAttitude;
}

void SetupRTIRunEvent()
{
	hRTIEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
}
void CloseRTIRunEvent(void)
{
	CloseHandle(hRTIEvent);
}

signed char WaitEvent(const unsigned char ucStepNo,const unsigned long ulMillisecond )
{
	signed char cReturn = 0;
	if( ucSynchronisedNo[ucStepNo] == 0 ) return 0;
	cReturn = ( WaitForSingleObject(hRTIEvent,ulMillisecond) == WAIT_OBJECT_0 );
	ResetEvent(hRTIEvent);
	return cReturn; 
}
void	SetRTIEvent(void)
{
	SetEvent(hRTIEvent);
}
void SendAgentRTIRegister(const unsigned char ucID,const unsigned char ucStepNo[MAX_STEP_NO])
{
	char chrSendMessage[MAX_RTI_MESSAGE_LENGTH];
	unsigned char ucLength = 0;
	unsigned char ucHead = 0;
	unsigned char i= 0;
	for(i=0;i<MAX_STEP_NO;i++)	ucSynchronisedNo[i] = ucStepNo[i];
	ucHead = ucIDRegisteProcessThread;
	ucRTIFederalNo = ucID;
	memcpy(chrSendMessage+ucLength,&ucHead,1);		ucLength += 1;
	memcpy(chrSendMessage+ucLength,&ucRTIFederalNo,1);		ucLength += 1;
	memcpy(chrSendMessage+ucLength,ucSynchronisedNo,MAX_STEP_NO);	ucLength += MAX_STEP_NO;
	chrSendMessage[ucLength] = 0;	SendRTIMessage(ulRTIChannelNo,chrSendMessage,ucLength);
}
void FeedbackAgentRTIReady(const unsigned char ucStepNo)
{
	char chrSendMessage[MAX_RTI_MESSAGE_LENGTH];
	unsigned char ucLength = 0;
	unsigned char ucHead = 0;
	if ( ucSynchronisedNo[ucStepNo] == 0 ) return;
	ucHead = ucIDStepFinished;
	memcpy(chrSendMessage+ucLength,&ucHead,1);		ucLength += 1;
	memcpy(chrSendMessage+ucLength,&ucSynchronisedNo[ucStepNo],1);		ucLength += 1;
	chrSendMessage[ucLength] = 0;		SendRTIMessage(ulRTIChannelNo,chrSendMessage,ucLength);
}
void SetQuitCommand(void);
unsigned short RTICallback(const char chrMessage[],const unsigned short usMessageLength)
{
	unsigned char ucHead = 0;
	unsigned char ucLength = 0;
	unsigned char ucID = 0;

	if ( usMessageLength < 2 ) return 0;
	memcpy(&ucHead,chrMessage,1);				ucLength += 1;
	switch (ucHead)
	{
	case ucIDRegisteThreadReturn:
		memcpy(&ucID,chrMessage+ucLength,1);	ucLength += 1;
		if ( ucID == ucRTIFederalNo )
			memcpy(ucSynchronisedNo,&chrMessage[ucLength],MAX_STEP_NO);	
		break;
	case ucIDStartNextStep:
		memcpy(&ucID,chrMessage+ucLength,1);	ucLength += 1;
		if ( ucSynchronisedNo[ucID] != 0 )
			SetEvent(hRTIEvent);
		break;
	case ucIDStopSimulation:
		SetQuitCommand();
		break;
	case ucIDIgnorAttitude:
		scIgnoreAttitude = 1;
		break;
	case  ucIDCalculateAtiitude:
		scIgnoreAttitude = 0;
		break;
	}
	return 0;
}