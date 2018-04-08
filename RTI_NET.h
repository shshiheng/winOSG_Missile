#ifndef __RTI_NET_H
#define __RTI_NET_H

#include "RTIConfig.h"

struct SRTIMessage
{
	unsigned short	usLength;
	char			chrMessage[MAX_RTI_MESSAGE_LENGTH];
};

unsigned short	CollectRTIMessages(const unsigned long ulChannelNo,struct SRTIMessage stcBufferOutput[]);
signed char		OpenRTIDevice(const unsigned char ucFederalNodeNo);
void			CloseRTIDevice(void);
signed char		SendRTIMessage(const unsigned long ulChannelNo, const char chrMessage[],const unsigned short usLength);
void			RegisteBUS(const unsigned long ulBUSNo,unsigned short (* pProcessFunction)(const char chrMessages[],const unsigned short usLength));

void AddToSendList(const unsigned long ulChannelNo,const char chrMessages[],const unsigned short usLength);

#endif
