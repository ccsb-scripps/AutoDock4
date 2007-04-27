/*

 $Id: readfield.h,v 1.5 2007/04/27 06:01:51 garrett Exp $

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

#ifndef READFIELD
#define READFIELD

#include "constants.h"
#include "openfile.h"
#include "stop.h"
#include "structs.h"

/*
void    readfield( Real *P_inv_spacing, 
                Real *P_spacing, 
                char  gdfldFileName[MAX_CHARS], 
                char  gpfFileName[MAX_CHARS], 
                int   gridpts1[SPACE], 
                int   gridpts[SPACE], 
		Real *xhi,
		Real *yhi,
		Real *zhi,
                Clock jobStart, 
                char  line[LINE_LEN], 
                Real *xlo, 
                Real *ylo, 
                Real *zlo, 
                char  macromolFileName[MAX_CHARS], 
                Real maP_center[SPACE], 
		struct tms tms_jobStart );
*/

void readfield( GridMapSetInfo *info, // *ptr_map_set_info
                char line[LINE_LEN],
                Clock jobStart,
                struct tms tms_jobStart );


#endif
