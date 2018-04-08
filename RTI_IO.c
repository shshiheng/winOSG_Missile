#include <stdio.h>
#include <string.h>

#include "RTI_NET.h"

#define USE_WINDOWS_PLATFORM
#ifdef USE_WINDOWS_PLATFORM
#include <windows.h>
static CRITICAL_SECTION gSection;
#endif

union UChannelNumber
{
	unsigned long ulNo;
	struct SChannelNumber
	{
		unsigned int usDeviceNo:16;
		unsigned int usChannelNo:16;
	}	stcNo;
};

#define MAX_PROGRAM_CHANNEL_NO (MAX_FEDERATE_NO*MAX_NET_NO+2*(MAX_CHANNEL_NO-MAX_NET_NO))

volatile static	struct SRTIMessage stcRTIBuffers[MAX_PROGRAM_CHANNEL_NO][RTI_BUFFER_SIZE];
volatile static unsigned short usRTIBufferStart[MAX_PROGRAM_CHANNEL_NO]={0}, usRTIBufferEnd[MAX_PROGRAM_CHANNEL_NO]={0};
unsigned short (*pFunctionArray[MAX_PROGRAM_CHANNEL_NO])(const char chrMessages[],const unsigned short usLength) = {NULL};
static unsigned short usTotalBusCount = 0;
static unsigned short usTotalNetCount = 0;

static unsigned short usRegistedDeviceNo[MAX_DEVICE_NO] = {0};
static unsigned short usTotalRegistedDeviceNo = 0;

static short sRegistedMessageBufferNo[MAX_CHANNEL_NO][MAX_DEVICE_NO] = {0};

void SetupRTIBuffer(void)
{
	unsigned short i = 0;
	unsigned short j = 0;
	for(i=0;i<MAX_PROGRAM_CHANNEL_NO;i++)
	{
		usRTIBufferStart[i] = 0;
		usRTIBufferEnd[i] = 0;
		pFunctionArray[i] = NULL;
	}

	for(i=0;i<MAX_DEVICE_NO;i++)	usRegistedDeviceNo[i] = 0;
	for(i=0;i<MAX_CHANNEL_NO;i++)	for(j=0;j<MAX_DEVICE_NO;j++)	sRegistedMessageBufferNo[i][j] = -1;

#ifdef USE_WINDOWS_PLATFORM
	InitializeCriticalSection(&gSection);
#endif
}
void SaveRTIMessageToBuffers(const unsigned short usNo, const char chrMessages[], const unsigned short usLength)
{
	stcRTIBuffers[usNo][usRTIBufferEnd[usNo]].usLength = usLength;
	memcpy((void*)(stcRTIBuffers[usNo][usRTIBufferEnd[usNo]].chrMessage),chrMessages,usLength);
	usRTIBufferEnd[usNo] = usRTIBufferEnd[usNo] + 1;
	if ( usRTIBufferEnd[usNo] >= RTI_BUFFER_SIZE )  usRTIBufferEnd[usNo] = 0;
}
unsigned short SaveRTIMessage(const unsigned short usBusNo, const unsigned short usDeviceNo,const char chrMessages[],const unsigned short usLength)
{
	unsigned short usIsSave = 0;
	unsigned short i = 0;
	unsigned short usSavedDeviceNo = 0;
	unsigned short usCount = 0;

#ifdef USE_WINDOWS_PLATFORM
	EnterCriticalSection(&gSection);
#endif

	for(i=0;i<usTotalRegistedDeviceNo;i++)
	{
		usSavedDeviceNo = usRegistedDeviceNo[i];
		if ( usSavedDeviceNo == usDeviceNo ) continue;
		if ( sRegistedMessageBufferNo[usBusNo][usSavedDeviceNo] == -1 )	continue;
		usCount++;
		if ( pFunctionArray[sRegistedMessageBufferNo[usBusNo][usSavedDeviceNo]] != NULL )
			usIsSave = pFunctionArray[sRegistedMessageBufferNo[usBusNo][usSavedDeviceNo]](chrMessages,usLength);
		else
			usIsSave = 1;
		if ( usIsSave == 1 )
			SaveRTIMessageToBuffers(sRegistedMessageBufferNo[usBusNo][usSavedDeviceNo],chrMessages,usLength);
		if ( usBusNo > MAX_NET_NO && usCount > 0 )	break;
	}

#ifdef USE_WINDOWS_PLATFORM
	LeaveCriticalSection(&gSection);
#endif

	return usCount;
}
signed char BroadcastRTIMessage(const unsigned short ulChannelNo, const char chrSendBuffer[],const unsigned short usLength);
signed char SendRTIMessage(const unsigned long ulChannelNo, const char chrSendBuffer[],const unsigned short usLength)
{
	union UChannelNumber stcChannelNo;
	unsigned short usFind = 0;
	stcChannelNo.ulNo = ulChannelNo;
	usFind = SaveRTIMessage(stcChannelNo.stcNo.usChannelNo,stcChannelNo.stcNo.usDeviceNo, chrSendBuffer,usLength);
	if ( stcChannelNo.stcNo.usChannelNo > MAX_NET_NO && usFind > 0 ) 
		return 0;
	else
		return ( BroadcastRTIMessage(stcChannelNo.stcNo.usChannelNo, chrSendBuffer,usLength) );
}

unsigned short CollectRTIMessages(const unsigned long ulChannelNo, struct SRTIMessage stcBufferOutput[])
{
	unsigned short usLength=0;
	unsigned short usEnd = 0;
	unsigned short usStart = 0; 
	unsigned short usSize = sizeof(struct SRTIMessage);
	unsigned short usBusNo = 0;
	union UChannelNumber stcChannelNo;
	stcChannelNo.ulNo = ulChannelNo;	usBusNo = stcChannelNo.stcNo.usChannelNo ;

	usBusNo = sRegistedMessageBufferNo[stcChannelNo.stcNo.usChannelNo][stcChannelNo.stcNo.usDeviceNo] ; 
	if ( usBusNo == -1 )	return 0;

	usEnd = usRTIBufferEnd[usBusNo];
	usStart = usRTIBufferStart[usBusNo]; 
	if (usEnd == usStart) 
		return(0);
	else if (usEnd > usStart)
	{
		memcpy((void*)stcBufferOutput,(void*)(&stcRTIBuffers[usBusNo][usStart]),(usEnd-usStart)*usSize);
		usLength = usEnd-usStart;
	}
	else
	{
		memcpy((void*)stcBufferOutput,(void*)(&stcRTIBuffers[usBusNo][usStart]),(RTI_BUFFER_SIZE-usStart)*usSize);
		if ( usEnd != 0 )	memcpy((void*)(stcBufferOutput+RTI_BUFFER_SIZE-usStart),(void*)(&stcRTIBuffers[usBusNo][0]),usEnd*usSize);
		usLength = RTI_BUFFER_SIZE-usStart+usEnd;
	}
	usRTIBufferStart[usBusNo] = usEnd;
	return usLength;
}

void RegisteBUS(const unsigned long ulChannelNo,unsigned short (* pProcessFunction)(const char chrMessages[],const unsigned short usLength))
{
	unsigned short i = 0;
	short sTestDevice = -1;
	union UChannelNumber stcChannelNo;
	stcChannelNo.ulNo = ulChannelNo;

	if(usTotalBusCount > MAX_CHANNEL_NO || stcChannelNo.stcNo.usChannelNo > MAX_CHANNEL_NO - 1)		
	{
#ifdef DEBUG_ERROR_PRINTF
		printf("Channel No:%d exceed range!\n",ulBUSNo);	
#endif
		return;
	}

	if ( sRegistedMessageBufferNo[stcChannelNo.stcNo.usChannelNo][stcChannelNo.stcNo.usDeviceNo] != -1 )
	{
#ifdef DEBUG_ERROR_PRINTF
		printf("Repeat Register channel %ld to use\n",i);
#endif
		return;
	}

	for(i=0;i<usTotalRegistedDeviceNo;i++)
	{
		if ( usRegistedDeviceNo[i] == stcChannelNo.stcNo.usDeviceNo )
		{
			sTestDevice = i;
			break;
		}
	}
	if ( sTestDevice < 0 )
	{
		usRegistedDeviceNo[usTotalRegistedDeviceNo] = stcChannelNo.stcNo.usDeviceNo;
		usTotalRegistedDeviceNo++;
	}

	sRegistedMessageBufferNo[stcChannelNo.stcNo.usChannelNo][stcChannelNo.stcNo.usDeviceNo] = usTotalBusCount;
	pFunctionArray[usTotalBusCount] = pProcessFunction;
	usTotalBusCount ++;	

#ifdef DEBUG_ERROR_PRINTF
	printf("Register channel %ld to use\n",ulTotalBusCount);
#endif

}