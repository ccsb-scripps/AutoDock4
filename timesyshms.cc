/* timesyshms.cc */

    #include <stdio.h>
    #include <sys/types.h>
    #include <sys/times.h>
    #include <time.h>
    #include <unistd.h>
    #include "timesyshms.h"


extern  FILE    *logFile;
extern	float	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( Clock  duration,
		 struct tms *start,
		 struct tms *end )

/*----------------------------------------------------------------------------*/

{
    int   h,
          m;
    float t,
	  T,
	  s;
    const float min = 60.,
                hrs = 3600.;
 

    (void)fprintf( logFile, "Real= " );
    t = (float)duration * idct;
    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2fs",       s );
        else
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( logFile, ",  CPU= " );
    t = (float)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2fs",       s );
        else
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( logFile, ",  System= " );
    t = (float)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;
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
