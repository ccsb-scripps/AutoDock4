/*
 $Id: getInitialState.h,v 1.7.6.1 2005/10/10 16:47:29 alther Exp $
*/

#ifndef GETINITIALSTATE
#define GETINITIALSTATE

#include "autocomm.h"
#include "constants.h"
#include "grid.h"
#include "structs.h"
#include "typedefs.h"


void getInitialState(FloatOrDouble *Addr_e0,
                     FloatOrDouble e0max,
	                  State *sInit,
	                  State *sMin,
	                  State *sLast,

                     Boole B_RandomTran0,
                     Boole B_RandomQuat0,
                     Boole B_RandomDihe0,

                     FloatOrDouble charge[MAX_ATOMS],
                     FloatOrDouble abs_charge[MAX_ATOMS],
                     FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                     FloatOrDouble q1q2[MAX_NONBONDS],
                     FloatOrDouble crd[MAX_ATOMS][SPACE],
                     FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
                     char  atomstuff[MAX_ATOMS][MAX_CHARS],
                     FloatOrDouble elec[MAX_ATOMS],
                     FloatOrDouble emap[MAX_ATOMS],

                     EnergyTables *ptr_ad_energy_tables,

                     Boole B_calcIntElec,
                     FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                     int   natom,
                     int   Nnb,
                     int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                     int   ntor,
                     int   tlist[MAX_TORS][MAX_ATOMS],
                     int   type[MAX_ATOMS],
                     FloatOrDouble vt[MAX_TORS][SPACE],
                     int   irun1,
                     int   outlev,
	                  int   MaxRetries,

                     FloatOrDouble torsFreeEnergy,

                     int   ligand_is_inhibitor,

                     int ignore_inter[MAX_ATOMS],

                     const Boole         B_include_1_4_interactions,
                     const FloatOrDouble scale_1_4,

                     const ParameterEntry parameterArray[MAX_MAPS],

                     const FloatOrDouble unbound_internal_FE,

                     GridMapSetInfo *info
                    );

#endif
