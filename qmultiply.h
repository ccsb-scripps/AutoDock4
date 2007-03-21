
#ifndef QMULTIPLY
#define QMULTIPLY
#include <stdio.h>
#include "constants.h"
#include "structs.h"

Quat uniformQuat( void );
Quat convertQuatToRot( Quat q );
Quat convertRotToQuat( Quat q );
Quat raaToQuat( const Real raa[3], Real angle );
Quat normQuat( Quat q );
Quat normRot( Quat q );
Quat conjugate( const Quat q );
Quat inverse( const Quat q );
Quat slerp( const Quat q1, const Quat q2, const double u );
Quat axisRadianToQuat( const Real ax, const Real ay, const Real az, const Real angle );
Quat axisDegreeToQuat( const Real ax, const Real ay, const Real az, const Real angle );
Quat quatComponentsToQuat( const Real qx, const Real qy, const Real qz, const Real qw );

void qmultiply( Quat *q, register const Quat *ql, register const Quat *qr );
void qconjmultiply( Quat *q, register const Quat *ql, register const Quat *qr );
void  mkUnitQuat( Quat *q );
void printQuat_q( FILE *fp, Quat q );
void printQuat_r( FILE *fp, Quat q );
void printQuat( FILE *fp, Quat q );
void unitQuat2rotation( Quat *q );
void print_q_reorient_message( FILE *logFile, Quat q_reorient );
void create_random_orientation( Quat *ptr_quat );
void assertQuatOK( const Quat q );
#endif
