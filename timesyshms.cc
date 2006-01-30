/*

 $Id: timesyshms.cc,v 1.6 2006/01/30 23:08:48 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* timesyshms.cc */

#include <stdio.h>
#include <sys/types.h>

#ifndef _WIN32
#include <sys/times.h>
#include <unistd.h>
#else
#include "times.h"
#endif
#include "timesyshms.h"

#include <time.h>


extern  FILE    *logFile;
extern	FloatOrDouble	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( Clock     duration,
                 struct tms  *start,
                 struct tms  *end)

/*----------------------------------------------------------------------------*/

{
    int   h, m;
    FloatOrDouble t, T, s;
    const FloatOrDouble min = 60., hrs = 3600.;

    (void)fprintf( logFile, "Real= " );
    t = (FloatOrDouble)duration * idct;
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2lfs",       (double)s );
        else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, (double)s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, (double)s );
    }

    (void)fprintf( logFile, ",  CPU= " );
    t =      (FloatOrDouble)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2lfs",       (double)s );
        else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, (double)s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, (double)s );
    }

    (void)fprintf( logFile, ",  System= " );
    t = (FloatOrDouble)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2lfs",       (double)s );
        else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, (double)s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, (double)s );
    }

    (void)fprintf( logFile, "\n" );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
