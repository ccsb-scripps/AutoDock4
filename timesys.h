

#ifndef TIMESYS
#define TIMESYS
#include <sys/types.h>

#ifdef _WIN32
#include "times.h"
#else
#include <sys/times.h>
#endif

#include <time.h>
#include "autocomm.h"
void  timesys( Clock  duration,
               struct tms *start,
               struct tms *end );
#endif
