/*

 $Id: timesys.cc,v 1.2 2003/02/26 01:47:47 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* timesys.cc */


    #include <stdio.h>
    #include <unistd.h>
    #include "timesys.h"

extern  FILE    *logFile;
extern	FloatOrDouble	idct;

/*----------------------------------------------------------------------------*/

void timesys( Clock  duration,
	      struct tms *start,
	      struct tms *end )

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (FloatOrDouble)duration * idct,
                         (FloatOrDouble)(end->tms_utime  - start->tms_utime) * idct,
                         (FloatOrDouble)(end->tms_stime  - start->tms_stime) * idct );
}
/* EOF */
