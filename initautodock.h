#ifndef INITAUTODOCK
#define INITAUTODOCK

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "print_2x.h"

void initautodock(
	    char  atomstuff[MAX_ATOMS][MAX_CHARS],
	    float crd[MAX_ATOMS][SPACE],
	    float crdpdb[MAX_ATOMS][SPACE],
	    float xhi,
	    float yhi,
	    float zhi,
	    float xlo,
	    float ylo,
	    float zlo,
	    int   natom,
	    int   ntor,
	    State *s0,
	    int   tlist[MAX_TORS][MAX_ATOMS],
	    float vt[MAX_TORS][SPACE],
	    int   outlev);
#endif
