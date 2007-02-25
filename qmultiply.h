
#ifndef QMULTIPLY
#define QMULTIPLY
#include <stdio.h>
#include "constants.h"
#include "structs.h"

void qmultiply( Quat *q,
                Quat *ql,
                Quat *qr );
void  mkUnitQuat( Quat *q );
void printQuat( FILE *fp, Quat q );
Quat uniformQuat( void );
void unitQuat2rotation( Quat *q );
Quat convertQuatToRot( Quat q );
Quat convertRotToQuat( Quat q );
Quat normQuat( Quat q );
Quat normRot( Quat q );
void print_q_reorient_message( FILE *logFile, Quat q_reorient );
void create_random_orientation( Quat *ptr_quat );
#endif
