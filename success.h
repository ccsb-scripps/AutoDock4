
#ifndef SUCCESS
#define SUCCESS
#include "constants.h"
#include "print_2x.h"
#include "timesyshms.h"
void  success( char  hostnm[MAX_CHARS],
               Clock jobStart,
               struct tms tms_jobStart );
#endif
