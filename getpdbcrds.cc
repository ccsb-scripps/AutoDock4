/*

 $Id: getpdbcrds.cc,v 1.3 2006/04/25 22:32:14 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* getpdbcrds.cc */

    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    #include <ctype.h>
    #include "getpdbcrds.h"


extern FILE *logFile;
extern char *programname;


int getpdbcrds( char rms_ref_crds_FN[MAX_CHARS],
		Real ref_crds[MAX_ATOMS][SPACE] )
{
    int ii=0;
    int natoms=0;
    char line[LINE_LEN];
    char str[4][WORDLEN];
    char rec5[5];
    FILE *rms_ref_FilePtr;

    if ( !openfile( rms_ref_crds_FN, "r", &rms_ref_FilePtr )) {
	fprintf( logFile, "%s: ERROR!  Sorry, could not open file \"%s\" for reading.\n", programname,  rms_ref_crds_FN );
	return -1;
    }

    pr (logFile, "\nRMS Reference Coordinates from \"%s\":-\n\n", rms_ref_crds_FN);

    while ( fgets(line, LINE_LEN, rms_ref_FilePtr) != NULL ) {

	for (ii = 0; ii < 4; ii++) {
	    rec5[ii] = (char)tolower( (int)line[ii] );
	}
	if (equal(rec5,"atom",4) || equal(rec5,"heta",4)) {

	    if (natoms < MAX_ATOMS) {
		sscanf( &line[30], "%s %s %s", str[X], str[Y], str[Z] );
		ref_crds[natoms][X] = atof( str[X] );
		ref_crds[natoms][Y] = atof( str[Y] );
		ref_crds[natoms][Z] = atof( str[Z] );
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
