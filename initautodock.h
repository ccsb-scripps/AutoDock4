#ifndef INITAUTODOCK
#define INITAUTODOCK

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "print_2x.h"

void initautodock(
	    char  atomstuff[MAX_ATOMS][MAX_CHARS],
	    FloatOrDouble crd[MAX_ATOMS][SPACE],
	    FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
	    int   natom,
	    int   ntor,
	    State *s0,
	    int   tlist[MAX_TORS][MAX_ATOMS],
	    FloatOrDouble vt[MAX_TORS][SPACE],
	    int   outlev,
            GridMapSetInfo *info );
#endif
