/* usage.cc */

    #include <stdio.h>
    #include "usage.h"


extern char    *programname;
extern char    AutoDockHelp[];

/*----------------------------------------------------------------------------*/
void usage( void )
/*----------------------------------------------------------------------------*/
{
    fprintf(stderr, "usage: %s %s\n", programname, AutoDockHelp);
}
/* EOF */
