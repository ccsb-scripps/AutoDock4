/*

 $Id: qmultiply.h,v 1.12 2009/05/08 23:02:16 rhuey Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

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
Real quatDifferenceToAngle( const Quat ql, const Quat qr );
Real quatDifferenceToAngleDeg( const Quat ql, const Quat qr );
Quat conjugate( const Quat q );
Quat inverse( const Quat q );
Quat slerp( const Quat qa, const Quat qb, const double t );
Quat slerp0( const Quat qa, const Quat qb, const double t );
Quat slerp1( const Quat qa, const Quat qb, const double t );
Quat axisRadianToQuat( const Real ax, const Real ay, const Real az, const Real angle );
Quat axisDegreeToQuat( const Real ax, const Real ay, const Real az, const Real angle );
Quat quatComponentsToQuat( const Real qx, const Real qy, const Real qz, const Real qw );

void qmultiply( Quat *q, register const Quat *ql, register const Quat *qr );
void qconjmultiply( Quat *q, register const Quat *ql, register const Quat *qr );
void mkUnitQuat( Quat *q );
void printQuat_q( FILE *fp, Quat q );
void printQuat_r( FILE *fp, Quat q );
void printQuat( FILE *fp, Quat q );
void debugQuat( FILE *fp, Quat q, unsigned int linenumber, char *message );
Quat uniformQuatByAmount( Real amount );
void unitQuat2rotation( Quat *q );
void print_q_reorient_message( FILE *logFile, Quat q_reorient );
void create_random_orientation( Quat *ptr_quat );
//void assertQuatOK( const Quat q );
const Quat identityQuat();
Real a_range_reduction( Real a );
Real alerp( Real a, Real b, Real fract );
#endif
