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
            FloatOrDouble *Addr_e0,
            FloatOrDouble e0max,

	    State *sInit,
	    State *sMin,
	    State *sLast,

            Boole B_RandomTran0,
            Boole B_RandomQuat0,
            Boole B_RandomDihe0,

            FloatOrDouble charge[MAX_ATOMS],
            FloatOrDouble q1q2[MAX_NONBONDS],
            FloatOrDouble crd[MAX_ATOMS][SPACE],
            FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
            char  atomstuff[MAX_ATOMS][MAX_CHARS],
            FloatOrDouble elec[MAX_ATOMS],
            FloatOrDouble emap[MAX_ATOMS],
            FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
            Boole B_calcIntElec,
            FloatOrDouble xhi,
            FloatOrDouble yhi,
            FloatOrDouble zhi,
            FloatOrDouble xlo,
            FloatOrDouble ylo,
            FloatOrDouble zlo,
            FloatOrDouble inv_spacing,
            FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int   natom,
            int   Nnb,
            int   nonbondlist[MAX_NONBONDS][4],
            int   ntor,
            int   tlist[MAX_TORS][MAX_ATOMS],
            int   type[MAX_ATOMS],
            FloatOrDouble vt[MAX_TORS][SPACE],
            int   irun1,
            int   outlev,
	       int   MaxRetries,
	       FloatOrDouble torsFreeEnergy,
            int   ligand_is_inhibitor,
            int ignore_inter[MAX_ATOMS]);
#endif
