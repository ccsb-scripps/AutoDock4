/* test_times.cc */

// possibly unnecessary // #include <iostream.h>
    #include <stdio.h>
#include <sys/types.h> // time_t time(time_t *tloc);
#include <time.h>      // time_t time(time_t *tloc);
#include <sys/times.h>
#include <unistd.h> // sysconf
#include <stdlib.h> // exit
#include "timesyshms.h"

#ifdef __alpha
#define Clock time_t
#else
#define Clock clock_t
#endif /* #ifdef __alpha */

float idct;


#ifdef USE_INT_AS_LONG
    typedef int  FourByteLong;
    typedef unsigned int UnsignedFourByteLong;
#else
    typedef long FourByteLong;
    typedef unsigned long UnsignedFourByteLong;
#endif


int main( int argc, char **argv, char **envp );

int main( int argc, char **argv, char **envp )
{
    static FourByteLong clktck = 0;
    struct tms tms_jobStart;
    struct tms tms_jobEnd;
    Clock  jobStart;
    Clock  jobEnd;
    long i=0L, j=0L;

    if (clktck == 0) {        /* fetch clock ticks per second first time */
        if ( (clktck = sysconf(_SC_CLK_TCK)) < (FourByteLong)0L) {
            (void) printf("\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
            exit( -1 );
        } else {
            idct = (float)1. / (float)clktck;
            (void) printf("\n\nFYI:  Number of clock ticks per second = %d\nFYI:  Elapsed time per clock tick = %.3e seconds\n\n\n\n", clktck, idct);
        }
    }

    jobStart = times( &tms_jobStart );

    for (i=0; i<1e8; i++) {
      j = i;
      }

/*
** Get the time at the start of the run...
*/
    jobEnd = times( &tms_jobEnd );
    (void) printf( "\nRun completed;  time taken for this run:\n");
    timesyshms( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd );


}

    #include <sys/types.h>
    #include <sys/times.h>
    #include <time.h>
    #include <unistd.h>
    #include "timesyshms.h"


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
 

    (void)fprintf( stdout, "Real= " );
    t = (float)duration * idct;
    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(stdout,       "%.2fs",       s );
        else
            (void)fprintf(stdout,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(stdout, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( stdout, ",  CPU= " );
    t = (float)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(stdout,       "%.2fs",       s );
        else
            (void)fprintf(stdout,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(stdout, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( stdout, ",  System= " );
    t = (float)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((float)h)*hrs;
    m = (int)(T/min);
    s = T - ((float)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(stdout,       "%.2fs",       s );
        else
            (void)fprintf(stdout,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(stdout, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( stdout, "\n" );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
