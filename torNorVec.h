
#ifndef TORNORVEC
#define TORNORVEC
#include "constants.h"
#include "stop.h"
void  torNorVec( Real crdpdb[MAX_ATOMS][SPACE],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 Real vt[MAX_TORS][SPACE] );
#endif
