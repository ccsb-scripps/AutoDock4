/*

 $Id: check_header_int.cc,v 1.3 2007/04/27 06:01:48 garrett Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "check_header_int.h"


extern char *programname;
extern FILE *logFile;


void check_header_int( int i1, 
		       int i2, 
		       char axis, 
		       char filename[MAX_CHARS] )

{
    char message[LINE_LEN];

    if ( i1 != i2 ) { 

	sprintf( message, "%s: Wrong number of %c grid-points in grid-map file \"%s\".\n", programname, (char)axis, filename );

        print_2x( logFile, stderr, message );

	sprintf( message, "%s: Use either %d or %d throughout!\n", programname, i1, i2 );

        print_2x( logFile, stderr, message );

        stop("Using wrong grid-map file.\n");
    }
}
/* EOF */
