/*

 $Id: check_header_line.cc,v 1.2 2003/02/26 00:37:28 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* check_header_line.cc */

    #include <stdio.h>
    #include <string.h>
    #include "check_header_line.h"

extern char *programname;
extern FILE *logFile;


void check_header_line( char s1[], char s2[] )

{
    if ( !equal(s1, s2, strlen(s1) ) ) {

	fprintf( logFile,"%s: Filename mismatch:\n\t\t\"%s\" :: \"%s\"\n", programname, s1,s2); 
    }
}
