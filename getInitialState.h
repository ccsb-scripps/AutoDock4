#ifndef GETINITIALSTATE
#define GETINITIALSTATE

#include "constants.h"
#include "qmultiply.h"
#include "stateLibrary.h"
#include "initautodock.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "cnv_state_to_coords.h"
#include "prInitialState.h"
#include "timesys.h"

void getInitialState(  
            float *Addr_e0,
            float e0max,

	    State *sInit,
	    State *sMin,
	    State *sLast,

            Boole B_RandomTran0,
            Boole B_RandomQuat0,
            Boole B_RandomDihe0,

            float charge[MAX_ATOMS],
            float q1q2[MAX_NONBONDS],
            float crd[MAX_ATOMS][SPACE],
            float crdpdb[MAX_ATOMS][SPACE],
            char  atomstuff[MAX_ATOMS][MAX_CHARS],
            float elec[MAX_ATOMS],
            float emap[MAX_ATOMS],
            float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
            Boole B_calcIntElec,
            float xhi,
            float yhi,
            float zhi,
            float xlo,
            float ylo,
            float zlo,
            float inv_spacing,
            float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int   natom,
            int   Nnb,
            int   nonbondlist[MAX_NONBONDS][2],
            int   ntor,
            int   tlist[MAX_TORS][MAX_ATOMS],
            int   type[MAX_ATOMS],
            float vt[MAX_TORS][SPACE],
            int   irun1,
            int   outlev,
	    int   MaxRetries,
	    float torsFreeEnergy,
      int   ligand_is_inhibitor);
#endif
