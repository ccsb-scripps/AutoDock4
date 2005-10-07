/*

 $Id: success.cc,v 1.2 2003/02/26 01:42:48 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* success.cc */

    #include <stdio.h>
    #include <sys/times.h>
    #include <sys/types.h>
    #include <time.h>
    #include "success.h"
    #include "timesyshms.h"


extern char *programname;
extern FILE *logFile;

void success( char hostnm[MAX_CHARS],
		Clock jobStart,
		struct tms tms_jobStart )

{
    char message[LINE_LEN];
    Clock jobEnd;
    struct tms tms_jobEnd;

    pr_2x( logFile, stderr, "\n" );
    pr_2x( logFile, stderr, UnderLine );
    prStr( message, "%s: Successful Completion on \"%s\"\n\n", programname, hostnm );
    pr_2x( logFile, stderr, message );

    jobEnd = times( &tms_jobEnd );

    timesyshms( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd );

    pr_2x( logFile, stderr, UnderLine );
}
/* EOF */
