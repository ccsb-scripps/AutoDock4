/*

 $Id: usage.cc,v 1.2 2003/02/26 01:50:51 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
