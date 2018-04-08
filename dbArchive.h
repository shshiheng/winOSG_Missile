#ifndef __DATABASE_ARCHIVE_H
#define __DATABASE_ARCHIVE_H

void AddDBDoubleData(const double DataValue);
void AddDBIntData(const long lDataValue);
void CreateDBMessage(const char chrName[]);
void SendDBMessage(void);

#endif