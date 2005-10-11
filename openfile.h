

#ifndef OPENFILE
#define OPENFILE

#include "constants.h"
#include "timesys.h"
#include "print_2x.h"
#include <sys/types.h>      /*time_t time(time_t *tloc); */
#include <time.h>           /*time_t time(time_t *tloc); */

#ifdef _WIN32
#include "times.h"
#else
#include <sys/times.h>
#endif

int  openfile( char  filename[MAX_CHARS],
               char  mode[],
               FILE  **fp );


int openFile( char       filename[MAX_CHARS],
              char       mode[],
              FILE       **fp,
              Clock      start,
              struct tms tms_start,
	      Boole	 mayExit);

FILE *ad_fopen(const char *path, const char *mode);

#endif
