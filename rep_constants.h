/*

 $Id: rep_constants.h,v 1.2 2007/04/27 06:01:51 garrett Exp $

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

#ifndef _REP_CONSTANTS
#define _REP_CONSTANTS

// Translation
#define X_TRANSLATION_INDEX 0
#define Z_TRANSLATION_INDEX 2

// Quaternion
#define QX_ROTATION_INDEX 3
#define QY_ROTATION_INDEX 4
#define QZ_ROTATION_INDEX 5
#define QW_ROTATION_INDEX 6

// Axis-Angle
#define X_ROTATION_INDEX 3
#define Y_ROTATION_INDEX 4
#define Z_ROTATION_INDEX 5

#define ROTATION_ANGLE_INDEX QW_ROTATION_INDEX

#define is_translation_index(i) (((i) >= X_TRANSLATION_INDEX) && ((i) <= Z_TRANSLATION_INDEX))
#define is_axis_index(i) (((i) >= X_ROTATION_INDEX) && ((i) <= Z_ROTATION_INDEX))
#define is_angle_index(i) ((i) == ROTATION_ANGLE_INDEX)
#define is_rotation_index(i) (((i) >= QX_ROTATION_INDEX) && ((i) <= QW_ROTATION_INDEX))

#endif
