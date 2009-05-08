/*

 $Id: check_header_line.cc,v 1.5 2009/05/08 23:02:11 rhuey Exp $

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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include "check_header_line.h"

extern char *programname;
extern FILE *logFile;


void check_header_line( char s1[], char s2[] )

{
    if ( !equal(s1, s2, strlen(s1) ) ) {

	    fprintf( logFile,"%s: Filename mismatch:\t\t\"%s\" :: \"%s\"\n", programname, s1,s2); 
    }
}
