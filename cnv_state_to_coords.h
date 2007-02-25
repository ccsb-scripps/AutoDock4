

#ifndef CNV_STATE_TO_COORDS
#define CNV_STATE_TO_COORDS

#include "constants.h"
#include "torsion.h"
#include "qtransform.h"

void cnv_state_to_coords( const State now,
                          Real vt[MAX_TORS][SPACE],
                          int tlist[MAX_TORS][MAX_ATOMS],
                          const int ntor,
                          Real crdpdb[MAX_ATOMS][SPACE],
                          Real crd[MAX_ATOMS][SPACE],
                          const int natom);
#endif
