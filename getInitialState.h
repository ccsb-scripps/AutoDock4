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
            Real *Addr_e0,
            Real e0max,

	    State *sInit,
	    State *sMin,
	    State *sLast,

            Boole B_RandomTran0,
            Boole B_RandomQuat0,
            Boole B_RandomDihe0,

            Real charge[MAX_ATOMS],
            Real abs_charge[MAX_ATOMS],
            Real qsp_abs_charge[MAX_ATOMS],
            Real q1q2[MAX_NONBONDS],
            Real crd[MAX_ATOMS][SPACE],
            Real crdpdb[MAX_ATOMS][SPACE],
            char  atomstuff[MAX_ATOMS][MAX_CHARS],
            Real elec[MAX_ATOMS],
            Real emap[MAX_ATOMS],

            EnergyTables *ptr_ad_energy_tables,

            Boole B_calcIntElec,
            Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int   natom,
            int   Nnb,
            int   **nonbondlist,
            int   ntor,
            int   tlist[MAX_TORS][MAX_ATOMS],
            int   type[MAX_ATOMS],
            Real vt[MAX_TORS][SPACE],
            int   irun1,
            int   outlev,
	        int   MaxRetries,

	        Real torsFreeEnergy,

            int   ligand_is_inhibitor,

            int ignore_inter[MAX_ATOMS],

            const Boole         B_include_1_4_interactions,
            const Real scale_1_4,

            const ParameterEntry parameterArray[MAX_MAPS],

            const Real unbound_internal_FE,

            GridMapSetInfo *info
           );

#endif
