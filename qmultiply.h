/*

 $Id: qmultiply.h,v 1.13 2010/08/27 00:05:08 mp Exp $

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
Quat convertQuatToRot( /* not const */ Quat q );
Quat convertRotToQuat( const Quat q );
Quat raaToQuat( const Real raa[3], Real angle );
Quat normQuat( /* not const */ Quat q );
Quat normRot( /* not const */ Quat q );
Real quatDifferenceToAngle( const Quat ql, const Quat qr );
Real quatDifferenceToAngleDeg( const Quat ql, const Quat qr );
Quat conjugate( const Quat q ) ;
Quat inverse( const Quat q ) ;
Quat slerp( const Quat qa, const Quat qb, const double t );
Quat slerp0( const Quat qa, const Quat qb, const double t );
Quat slerp1( const Quat qa, const Quat qb, const double t );
Quat axisRadianToQuat( const Real ax, const Real ay, const Real az, const Real angle );
Quat axisDegreeToQuat( const Real ax, const Real ay, const Real az, const Real angle );
Quat quatComponentsToQuat( const Real qx, const Real qy, const Real qz, const Real qw );

void qmultiply( Quat *const q, register const Quat *const ql, register const Quat *const qr );
void qconjmultiply( Quat *const q, register const Quat *const ql, register const Quat *const  qr );
void mkUnitQuat( Quat *const q );
void printQuat_q( FILE *const fp, const Quat q );
void printQuat_r( FILE *const fp, const Quat q );
void printQuat( FILE *const fp, const Quat q );
void debugQuat( FILE *const fp, const Quat q, const unsigned int linenumber, const char *const message );
Quat uniformQuatByAmount( const Real amount );
void unitQuat2rotation( /* not const */ Quat *q );
void print_q_reorient_message( FILE *const logFile, const Quat q_reorient );
void create_random_orientation( /* not const */ Quat *const ptr_quat );
//void assertQuatOK( const Quat q );
const Quat identityQuat() ;
Real a_range_reduction( /* not const */ Real a ) ;
Real alerp( /* not const */ Real a, /* not const */ Real b, const Real fract ) ;
#endif
