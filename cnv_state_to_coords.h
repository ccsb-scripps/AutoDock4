

#ifndef CNV_STATE_TO_COORDS
#define CNV_STATE_TO_COORDS

#include "constants.h"
#include "torsion.h"
#include "qtransform.h"

void cnv_state_to_coords( State now,
			  float vt[MAX_TORS][SPACE],
			  int   tlist[MAX_TORS][MAX_ATOMS],
			  int ntor,
			  float crdpdb[MAX_ATOMS][SPACE],
			  float crd[MAX_ATOMS][SPACE],
			  int   natom);
#endif
