/*

 $Id: print_rem.cc,v 1.4 2007/04/27 06:01:50 garrett Exp $

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
#include "print_rem.h"


void print_rem( FILE *outFile,
		int Rank,
		int NumMem,
		int Run,
		Real ref_rms)
{
    fprintf( outFile, "MODEL     %4d\n", Run );
    fprintf( outFile, "USER    Run = %d\n", Run );
    fprintf( outFile, "USER    Cluster Rank = %d\n", Rank );
    fprintf( outFile, "USER    Number of conformations in this cluster = %d\n", NumMem );
    fprintf( outFile, "USER  \n");
    fprintf( outFile, "USER    RMSD from reference structure       = %.3f A\n", ref_rms );
    fprintf( outFile, "USER  \n");
}
/* EOF */
