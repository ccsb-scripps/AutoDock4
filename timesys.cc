/* timesys.cc */


    #include <stdio.h>
    #include <unistd.h>
    #include "timesys.h"

extern  FILE    *logFile;
extern	float	idct;

/*----------------------------------------------------------------------------*/

void timesys( Clock  duration,
	      struct tms *start,
	      struct tms *end )

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (float)duration * idct,
                         (float)(end->tms_utime  - start->tms_utime) * idct,
                         (float)(end->tms_stime  - start->tms_stime) * idct );
}
/* EOF */
