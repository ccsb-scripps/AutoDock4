
#ifndef PRINT_AVSFLD
#define PRINT_AVSFLD
#include "constants.h"
void  print_avsfld(FILE  *logFile,
                   int   veclen,
                   int   natom,
                   int   nframe,
                   int   offset[VECLENMAX],
                   int   stride,
                   char  label[MAX_CHARS],
                   char  filename[MAX_CHARS] );
#endif
