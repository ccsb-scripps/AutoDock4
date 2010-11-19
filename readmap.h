/*

 $Id: readmap.h,v 1.7 2007/04/27 06:01:51 garrett Exp $

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


#ifndef MAPC2F
#define MAPC2F
#include "constants.h"
#include "openfile.h"
#include "warn_bad_file.h"
#include "strindex.h"
#include "print_2x.h"
#include "check_header_line.h"
#include "warn_bad_file.h"
#include "check_header_float.h"
#include "check_header_int.h"
#include "timesys.h"
Real   mapc2f( char C_mapValue );
#endif

#ifndef READMAP
#define READMAP
#include "constants.h"
#include "openfile.h"
#include "warn_bad_file.h"
#include "strindex.h"
#include "print_2x.h"
#include "check_header_line.h"
#include "warn_bad_file.h"
#include "check_header_float.h"
#include "check_header_int.h"
#include "timesys.h"
#include "structs.h"

Statistics readmap( char line[LINE_LEN],
             int outlev,

             Clock jobStart,
             struct tms tmsJobStart,
        
             Boole B_charMap,

             Boole *P_B_HaveMap, 
             int num_maps, 
             
             GridMapSetInfo *info,
             Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
             // double *maps 
             char map_type
             );

#endif
