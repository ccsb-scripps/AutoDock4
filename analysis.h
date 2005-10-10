/*
 $Id: analysis.h,v 1.7.6.1 2005/10/10 16:41:40 alther Exp $
*/

#ifndef ANALYSIS
#define ANALYSIS

#include "autocomm.h"
#include "constants.h"
#include "grid.h"
#include "structs.h"
#include "typedefs.h"


void  analysis( int   Nnb,
                char  atomstuff[MAX_ATOMS][MAX_CHARS],
                FloatOrDouble charge[MAX_ATOMS],
                FloatOrDouble abs_charge[MAX_ATOMS],
                FloatOrDouble qsp_abs_charge[MAX_ATOMS],
	             Boole B_calcIntElec,
	             FloatOrDouble q1q2[MAX_NONBONDS],
                FloatOrDouble clus_rms_tol,
                FloatOrDouble crdpdb[MAX_ATOMS][SPACE],

                EnergyTables *ptr_ad_energy_tables,

                FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                FloatOrDouble econf[MAX_RUNS],
                int   irunmax,
                int   natom,
                int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                int   nconf,
                int   ntor,
                State hist[MAX_RUNS],
                char  smFileName[MAX_CHARS],
                FloatOrDouble sml_center[SPACE],
                Boole B_symmetry_flag,
                int   tlist[MAX_TORS][MAX_ATOMS],
                int   type[MAX_ATOMS],
                FloatOrDouble vt[MAX_TORS][SPACE],
		          char  rms_ref_crds[MAX_CHARS],
		          FloatOrDouble torsFreeEnergy,
                Boole B_write_all_clusmem,
                int ligand_is_inhibitor,
                Boole B_template,
                FloatOrDouble template_energy[MAX_ATOMS], // template energy value for each atom
                FloatOrDouble template_stddev[MAX_ATOMS], // and standard deviation of this energy
                int   outlev,
                int   ignore_inter[MAX_ATOMS],
                const Boole   B_include_1_4_interactions,
                const FloatOrDouble scale_1_4,

                const ParameterEntry parameterArray[MAX_MAPS],
                const FloatOrDouble unbound_internal_FE,

                GridMapSetInfo *info );
#endif
