/*

 $Id: timesyshms.cc,v 1.2 2003/02/26 01:47:59 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* timesyshms.cc */

    #include <stdio.h>
    #include <sys/types.h>
    #include <sys/times.h>
    #include <time.h>
    #include <unistd.h>
    #include "timesyshms.h"


extern  FILE    *logFile;
extern	FloatOrDouble	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( Clock  duration,
		 struct tms *start,
		 struct tms *end )

/*----------------------------------------------------------------------------*/

{
    int   h,
          m;
    FloatOrDouble t,
	  T,
	  s;
    const FloatOrDouble min = 60.,
                hrs = 3600.;
 

    (void)fprintf( logFile, "Real= " );
    t = (FloatOrDouble)duration * idct;
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2fs",       s );
        else
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( logFile, ",  CPU= " );
    t = (FloatOrDouble)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2fs",       s );
        else
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( logFile, ",  System= " );
    t = (FloatOrDouble)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2fs",       s );
        else
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( logFile, "\n" );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
