/*

 $Id: qmultiply.h,v 1.7 2007/04/27 06:01:50 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

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
Quat conjugate( const Quat q );
Quat inverse( const Quat q );
Quat slerp( const Quat q1, const Quat q2, const double u );
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
void assertQuatOK( const Quat q );
const Quat identityQuat();
#endif
