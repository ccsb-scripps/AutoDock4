
#ifndef QMULTIPLY
#define QMULTIPLY
#include "constants.h"
#include "structs.h"

void qmultiply( Quat *q,
		Quat *ql,
		Quat *qr );

void  mkUnitQuat( Quat *q );

void printQuat( Quat q );

#endif
