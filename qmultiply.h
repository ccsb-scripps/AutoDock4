/*

 $Id: qmultiply.h,v 1.15 2011/03/08 04:18:37 mp Exp $

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
Quat convertQuatToRot( const Quat& q );
Quat convertRotToQuat( const Quat& q );
Quat raaToQuat( const Real raa[3], ConstReal angle );
Quat normQuat( /* not const */ Quat q );
Quat normRot( /* not const */ Quat q );
Real quatDifferenceToAngle( const Quat& ql, const Quat& qr );
Real quatDifferenceToAngleDeg( const Quat& ql, const Quat& qr );
Quat conjugate( const Quat& q );
Quat inverse( const Quat& q );
Quat slerp( const Quat& qa, const Quat& qb, ConstDouble t );
Quat slerp0( const Quat& qa, const Quat& qb, ConstDouble t );
Quat slerp1( const Quat& qa, const Quat& qb, ConstDouble t );
Quat axisRadianToQuat( ConstReal ax, ConstReal ay, ConstReal az, ConstReal angle );
Quat axisDegreeToQuat( ConstReal ax, ConstReal ay, ConstReal az, ConstReal angle );
Quat quatComponentsToQuat( ConstReal qx, ConstReal qy, ConstReal qz, ConstReal qw );

void qmultiply( Quat *const q, register const Quat *const ql, register const Quat *const qr );
void qconjmultiply( Quat *const q, register const Quat *const ql, register const Quat *const qr );
void mkUnitQuat( Quat *const q );
void printQuat_q( FILE *const fp, const Quat& q );
void printQuat_r( FILE *const fp, const Quat& q );
void printQuat( FILE *const fp, const Quat& q );
void debugQuat( FILE *const fp, const Quat& q, const unsigned int linenumber, const char *const message );
Quat uniformQuatByAmount( ConstReal amount );
void unitQuat2rotation( /* not const */ Quat *q );
void print_q_reorient_message( FILE *const logFile, const Quat& q_reorient );
void create_random_orientation( /* not const */ Quat *const ptr_quat );
//void assertQuatOK( const Quat q );
const Quat identityQuat();
Real a_range_reduction( ConstReal a );
Real alerp( ConstReal a, ConstReal b, ConstReal fract );
#endif
