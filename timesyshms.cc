/*

 $Id: timesyshms.cc,v 1.8 2007/04/27 06:01:51 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
extern	Real	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( Clock     duration,
                 struct tms  *start,
                 struct tms  *end)

/*----------------------------------------------------------------------------*/

{
    int   h, m;
    Real t, T, s;
    const Real min = 60., hrs = 3600.;

    (void)fprintf( logFile, "Real= " );
    t = (Real)duration * idct;
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2lfs",       (double)s );
        else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, (double)s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, (double)s );
    }

    (void)fprintf( logFile, ",  CPU= " );
    t =      (Real)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(logFile,       "%.2lfs",       (double)s );
        else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, (double)s );
    } else {
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, (double)s );
    }

    (void)fprintf( logFile, ",  System= " );
    t = (Real)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
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
