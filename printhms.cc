/* printhms.cc */

    #include <stdio.h>
    #include "printhms.h"

extern FILE *logFile;

void printhms( float t )

{
    int   h,
          m;
    float T, s;
    float min = 60.,
	  hrs = 3600.;

    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;

    if (h == 0) {
        if (m == 0)
            fprintf(logFile,       "%.2fs",       s );
        else
            fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }
}
