/*

 $Id: getpdbcrds.cc,v 1.10 2018/07/31 23:26:51 mp Exp $

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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "getpdbcrds.h"

/* define length of the x,y,z fields and staring column (1-origin) */
#define PDBCOORDLEN 8  
#define PDBCOORDORG 31


extern char *programname;


int getpdbcrds( const char *const rms_ref_crds_FN,
		/* not const */ Real ref_crds[MAX_ATOMS][SPACE], FILE *logFile)
{
    int natoms=0;
    char line[LINE_LEN];
    char rectype[6+1];
    char str[3][PDBCOORDLEN+1];
    FILE *rms_ref_FilePtr;

    if ( !openfile( rms_ref_crds_FN, "r", &rms_ref_FilePtr, logFile)) {
	fprintf( logFile, "%s: ERROR!  Sorry, could not open file \"%s\" for reading.\n", programname,  rms_ref_crds_FN );
	return -1;
    }

    pr (logFile, "\nRMS Reference Coordinates from \"%s\":-\n\n", rms_ref_crds_FN);

    while ( fgets(line, LINE_LEN, rms_ref_FilePtr) != NULL ) {

	for (int ii = 0; ii < 6; ii++) {
	    rectype[ii] = (char)tolower( (int)line[ii] );
	}
	rectype[6]='\0';
	if (equal(rectype,"atom  ") || equal(rectype,"hetatm")) {

	    if (natoms < MAX_ATOMS) {
		// WRONG sscanf( &line[30], "%s %s %s", str[X], str[Y], str[Z] );
		for(int c=0;c<3;c++)  {
		  (void) strncpy( str[c], &line[PDBCOORDORG+c*PDBCOORDLEN-1], (size_t)PDBCOORDLEN );
		  str[c][PDBCOORDLEN]='\0';
		ref_crds[natoms][c] = atof( str[c] );
		}
		pr (logFile, "Atom %5d,  x,y,z = %8.3f %8.3f %8.3f\n", natoms+1, ref_crds[natoms][X], ref_crds[natoms][Y], ref_crds[natoms][Z]);
	    } else {
		fprintf( logFile, "%s: ERROR!  Sorry, too many atoms in file \"%s\"\n", programname,  rms_ref_crds_FN );
		return -1;
	    }
	    ++natoms;
	} /* End if "atom" or "heta" */

    } /* End while there's a line to read... */

    /* Close the reference coordinates file... */
    (void) fclose(rms_ref_FilePtr);

    return natoms;
}
/* EOF */
