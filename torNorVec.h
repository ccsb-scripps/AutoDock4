
#ifndef TORNORVEC
#define TORNORVEC
#include "constants.h"
#include "stop.h"
#include "structs.h"

void  torNorVec( Real crdpdb[MAX_ATOMS][SPACE],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 Real vt[MAX_TORS][SPACE] );

void update_torsion_vectors( Real crdpdb[MAX_ATOMS][SPACE],
                             int ntor,
                             int  tlist[MAX_TORS][MAX_ATOMS],
                             Real vt[MAX_TORS][SPACE],
                             Molecule *ligand,
                             int debug );

#endif
