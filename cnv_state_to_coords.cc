/* cnv_state_to_coords.cc */

#include <math.h>
#include <string.h>

#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

extern FILE *logFile;

void cnv_state_to_coords( State now,
                          float vt[MAX_TORS][SPACE],
                          int   tlist[MAX_TORS][MAX_ATOMS],
                          int   ntor,
                          float crdpdb[MAX_ATOMS][SPACE],
                          float crd[MAX_ATOMS][SPACE],
                          int   natom)

{
    // register int i;
    //register int XYZ;

#ifdef DEBUG
    if (checkState(&now)) {
        // only if the state is valid can we proceed...
#endif

    //  Setting back to the original PDB
    //  coordinates ensures that cumulative
    //  rounding errors do not occur.
    //  this memcpy call...
    (void)memcpy(crd, crdpdb, natom*3*sizeof(float));

    //  is about 100x faster than these nested for-loops...
    //for (i = 0;  i < natom;  i++) {
        //for (XYZ = 0;  XYZ < SPACE;  XYZ++) {
            //crd[i][XYZ] = crdpdb[i][XYZ];
        //} //XYZ
    //} //i

    //  Apply torsions, quaternion rigid-body rotation
    //  and translation...
    if (ntor > 0) {
      torsion( now, crd, vt, tlist, ntor );
    }
    qtransform( now.T, now.Q, crd, natom );

#ifdef DEBUG
    } else {
        (void)fprintf(logFile,"WARNING: About to convert state to coordinates, but there is a problem with state! (The function 'checkState' failed).\n");
        // exit(-1);
    } // endif checkState
#endif
}

/* EOF */
