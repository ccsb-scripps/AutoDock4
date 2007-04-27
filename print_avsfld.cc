/*

 $Id: print_avsfld.cc,v 1.4 2007/04/27 06:01:50 garrett Exp $

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
#include "print_avsfld.h"
#include "autocomm.h"


void print_avsfld( FILE *logFile,
		   int veclen,
		   int natom,
		   int nframe,
		   int offset[VECLENMAX],
		   int stride,
		   char label[MAX_CHARS],
		   char filename[MAX_CHARS] )
{
    int i;

    fprintf( logFile, "AVSFLD: # AVS field file\n" );
    fprintf( logFile, "AVSFLD: #\n" );
    fprintf( logFile, "AVSFLD: # Created by AutoDock\n" );
    fprintf( logFile, "AVSFLD: #\n" );
    fprintf( logFile, "AVSFLD: ndim=2           # number of dimensions in the field\n" );
    fprintf( logFile, "AVSFLD: nspace=1         # number of physical coordinates\n" );
    fprintf( logFile, "AVSFLD: veclen=%-9d # vector size\n", veclen );
    fprintf( logFile, "AVSFLD: dim1=%-11d # atoms\n", natom );
    fprintf( logFile, "AVSFLD: dim2=%-11d # conformations\n", nframe );
    fprintf( logFile, "AVSFLD: data=Real       # data type (byte,integer,Real,double)\n" );
    fprintf( logFile, "AVSFLD: field=uniform    # field coordinate layout\n" );
    fprintf( logFile, "AVSFLD: label= %s\n", label );
    for (i=0; i<veclen; i++) {
        fprintf( logFile, "AVSFLD: variable %d file = %s filetype = ascii offset = %d stride = %d\n", 
					   i+1,    filename,                   offset[i],    stride );
    }
    fprintf( logFile, "AVSFLD: # end of file\n\n" );
    fflush( logFile );
}
/* EOF */
