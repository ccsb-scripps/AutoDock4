
#ifndef TORNORVEC
#define TORNORVEC
#include "constants.h"
#include "stop.h"
void  torNorVec( float crdpdb[MAX_ATOMS][SPACE],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 float vt[MAX_TORS][SPACE] );
#endif
