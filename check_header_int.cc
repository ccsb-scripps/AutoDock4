/*

 $Id: check_header_int.cc,v 1.2 2003/02/26 00:36:47 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* check_header_int.cc */

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
