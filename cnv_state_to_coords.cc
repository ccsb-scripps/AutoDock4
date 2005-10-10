/*

 $Id: cnv_state_to_coords.cc,v 1.4.6.1 2005/10/10 16:44:36 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* cnv_state_to_coords.cc */

#include <string.h>

#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

#ifdef LSQFIT
#include "rotlsqfit.h"
#endif

extern FILE *logFile;
extern int true_ligand_atoms;

void cnv_state_to_coords( State now,
                          FloatOrDouble vt[MAX_TORS][SPACE],
                          int   tlist[MAX_TORS][MAX_ATOMS],
                          int   ntor,
                          FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
                          FloatOrDouble crd[MAX_ATOMS][SPACE],
                          int   natom)

{
#ifdef DEBUG
    if (checkState(&now)) { // only if the state is valid can we proceed...
#endif

    //  Setting back to the original PDB
    //  coordinates ensures that cumulative
    //  rounding errors do not occur.

    (void) memcpy( crd, crdpdb,  natom * 3 * sizeof(FloatOrDouble));

    //  memcpy is about 100x faster than these nested for-loops...
    //  for (i = 0;  i < natom;  i++) for (XYZ = 0;  XYZ < SPACE;  XYZ++) crd[i][XYZ] = crdpdb[i][XYZ];

    //  Apply torsions, if any
    if (ntor > 0) {
        torsion( now, crd, vt, tlist, ntor );
    }

#ifdef LSQFIT
    //  Optionally, least-squares fit the current conformation on the original input conformation
    static double crdfit[MAX_ATOMS][SPACE];
    // double weights[MAX_ATOMS];
    // double vector[3];
    // double matrix[3][3];
    // rotlsqfit(crd, crdpdb, crdfit, weights, natom, matrix, vector)
    (void) rotlsqfit( crd, crdpdb, crdfit, NULL, natom, NULL, NULL);
    (void) memcpy( crd, crdfit,  natom * 3 * sizeof(FloatOrDouble));
#endif

    //  Apply quaternion rigid-body rotation and translation...
    qtransform( now.T, now.Q, crd, true_ligand_atoms );

#ifdef DEBUG
    } else {
        (void)fprintf(logFile,"WARNING: About to convert state to coordinates, but there is a problem with state! (The function 'checkState' failed).\n");
        // exit(-1);
    } // endif checkState
#endif
}

/* EOF */
