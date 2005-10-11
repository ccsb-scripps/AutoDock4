/*
 $Id: timesyshms.h,v 1.1.1.1.8.1 2005/10/11 00:09:23 alther Exp $
*/

#ifndef TIMESYSHMS
#define TIMESYSHMS

#ifdef _WIN32
   #include "times.h"
#else
   #include <sys/times.h>
#endif

#include "autocomm.h"
void  timesyshms( Clock  duration,
                  struct tms *start,
                  struct tms *end );
#endif
