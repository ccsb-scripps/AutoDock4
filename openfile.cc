/*

 $Id: openfile.cc,v 1.2 2003/02/26 01:22:37 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* openfile.cc */

    #include <stdio.h>
    #include <stdlib.h>
    #include "openfile.h"


extern char *programname;
extern FILE *logFile;

/*----------------------------------------------------------------------------*/
int openfile( char filename[MAX_CHARS],
	      char mode[],
	      FILE **fp )

{
	if ( (*fp = fopen(filename, mode)) == NULL ) {
		fprintf(stderr, "\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
		fprintf(logFile,"\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
		return( FALSE );
	}
	return( TRUE );
}

/*----------------------------------------------------------------------------*/
int openFile( char       filename[MAX_CHARS],
	      char       mode[],
	      FILE       **fp,
	      Clock      start,
	      struct tms tms_start,
	      Boole      mayExit)

{
    Clock  jobEnd;
    struct tms tms_jobEnd;

    if ( (*fp = fopen(filename, mode)) == NULL ) {
	fprintf(stderr, "\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
	fprintf(logFile,"\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);

	jobEnd = times( &tms_jobEnd );
	timesys( jobEnd - start, &tms_start, &tms_jobEnd );
	pr_2x( logFile, stderr, UnderLine );

	if (mayExit) {
	    exit(-1); /* END PROGRAM */
	} else {
	    return( FALSE );
	}
    }
    return( TRUE );
}
/*----------------------------------------------------------------------------*/
