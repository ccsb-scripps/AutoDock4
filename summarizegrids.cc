/*

 $Id: summarizegrids.cc,v 1.4 2007/04/27 06:01:51 garrett Exp $

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
#include "summarizegrids.h"


extern FILE *logFile;

void summarizegrids( char  atm_typ_str[ATOM_MAPS],
		     Real mapmax[MAX_MAPS],
		     Real mapmin[MAX_MAPS],
		     int   num_all_maps,
		     int   num_atm_maps )
{
    /*
     Summarize in the log-file what was found:
    */

    int i = 0;

    fprintf( logFile, "\n" );
    fprintf( logFile, UnderLine );
    fprintf( logFile, "SUMMARY OF GRID MAPS' ENERGIES:\n_______________________________\n\n\t\t\tvan der Waals or H-bond\n\tGrid\tAtom\t    Potential Energy\n\tMap\tType\tMinimum\t\tMaximum\n\t___\t____\t_________\t_________\n" );

    for ( i=0; i < num_atm_maps ; i++ ) {

	fprintf( logFile, "\t %d\t %c\t", (i+1), atm_typ_str[i]);

	/*fprintf( logFile, "<DEBUG 2>\n");  */

	fprintf( logFile, "%6.2f\t\t", mapmin[i]); 

	fflush( logFile );

	/* fprintf( logFile, "<DEBUG 3>\n");  */

	fflush( logFile );

	fprintf( logFile, "%6.2e\n", mapmax[i]);

	fflush( logFile );

	/* fprintf( logFile, "<DEBUG 4>\n");  */

	fflush( logFile );
    }

    fprintf( logFile, "\n\t\t\t    Electrostatic\n\tGrid\tGrid\t    Potential Energy:\n\tMap\tType\tMinimum\t\tMaximum\n\t___\t____\t_________\t_________\n" );

    fprintf( logFile, "\t %d\telec\t%+8.2f\t%+8.2f\n\n", num_all_maps, mapmin[num_atm_maps],mapmax[num_atm_maps]);

    /* fprintf( logFile, UnderLine ); */

    fflush( logFile );
}
/* EOF */
