/*

 $Id: warn_bad_file.cc,v 1.2 2003/02/26 01:50:58 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* warn_bad_file.cc */


    #include <stdio.h>
    #include "warn_bad_file.h"


extern char *programname;
extern FILE *logFile;

/*----------------------------------------------------------------------------*/

void warn_bad_file( char filename[MAX_CHARS],
		    char message[LINE_LEN] )

/*----------------------------------------------------------------------------*/
{
    char out[LINE_LEN];

    sprintf( out, "%s: WARNING: possibly unsuitable/old version grid-map file %s\n", programname, filename );
    print_2x( logFile, stderr, out );
    sprintf( out, "%s: %s\n", programname, message );
    print_2x( logFile, stderr, out );
}
/* EOF */
