/* stop.cc */

    #include <stdio.h>
    #include <stdlib.h>
    #include "stop.h"

extern char *programname;
extern FILE *logFile;


/*----------------------------------------------------------------------------*/
void stop( char reason[LINE_LEN] )
/*----------------------------------------------------------------------------*/
{
    char message[LINE_LEN];

    if (logFile == stdout) {
	fprintf( logFile, "%s: ERROR: %s\n", programname, reason);
	fprintf( logFile, "%s: Aborting...\n\n", programname);
	fprintf( logFile, "%s: Unsuccessful Completion.\n\n", programname);
	fflush(  logFile  );
    } else {
	sprintf( message, "%s: ERROR: %s\n", programname, reason);
	print_2x( logFile, stderr, message );
	sprintf( message, "%s: Aborting...\n\n", programname);
	print_2x( logFile, stderr, message );
	sprintf( message, "%s: Unsuccessful Completion.\n\n", programname);
	print_2x( logFile, stderr, message );
	fflush(logFile);
	fflush(stderr);
    }

    exit(-1);
}
/* EOF */
