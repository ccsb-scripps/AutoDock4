/* check_header_float.cc */

    #include <stdio.h>
    #include "stop.h"
    #include "check_header_float.h"


extern char *programname;
extern FILE *logFile;


void check_header_float( float f1, float f2, char keyword[], char filename[] )

{
    if ( f1 != f2 ) { 
        fprintf(logFile, "Wrong %s in grid-map file \"%s\".\n", keyword, filename);
        fprintf(stderr, "%s: Wrong %s in grid-map file \"%s\".\n", programname, keyword, filename);

        fprintf(logFile, "Use either %.3f or %.3f throughout!\n\n", f1,f2);
        fprintf(stderr, "%s: Use either %.3f or %.3f throughout!\n", programname, f1,f2);

        stop("Using wrong grid-map file.\n");
    }
}
/* EOF */
