
#ifndef TORNORVEC
#define TORNORVEC
#include "constants.h"
#include "stop.h"
void  torNorVec( FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 FloatOrDouble vt[MAX_TORS][SPACE] );
#endif
