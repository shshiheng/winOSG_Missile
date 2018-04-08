#ifndef   RTI_RUNEVENT_H
#define   RTI_RUNEVENT_H

void  SetupRTIRunEvent(void);
void  CloseRTIRunEvent(void);
signed char		WaitEvent(const unsigned char ucStepNo,const unsigned long ulMillisecond);
void			SendAgentRTIRegister(const unsigned char ucID,const unsigned char ucStepNo[]);
void			FeedbackAgentRTIReady(const unsigned char ucStepNo);
unsigned short	RTICallback(const char chrMessage[], const unsigned short usMessageLength);
void			SetRTIEvent(void);
#endif

