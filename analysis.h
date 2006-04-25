#include "constants.h"
#include "getpdbcrds.h"
#include "stateLibrary.h"
#include "cnv_state_to_coords.h"
#include "sort_enrg.h"
#include "cluster_analysis.h"
#include "prClusterHist.h"
#include "getrms.h"
#include "eintcal.h"
#include "trilinterp.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"

#ifndef ANALYSIS
#define ANALYSIS

void  analysis( int   Nnb, 
                char  atomstuff[MAX_ATOMS][MAX_CHARS], 
                Real charge[MAX_ATOMS], 
                Real abs_charge[MAX_ATOMS], 
                Real qsp_abs_charge[MAX_ATOMS], 
	            Boole B_calcIntElec,
	            Real q1q2[MAX_NONBONDS],
                Real clus_rms_tol, 
                Real crdpdb[MAX_ATOMS][SPACE], 
                
                EnergyTables *ptr_ad_energy_tables,

                Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                Real econf[MAX_RUNS], 
                int   irunmax, 
                int   natom, 
                int   **nonbondlist, 
                int   nconf, 
                int   ntor, 
                State hist[MAX_RUNS], 
                char  smFileName[MAX_CHARS], 
                Real sml_center[SPACE], 
                Boole B_symmetry_flag, 
                int   tlist[MAX_TORS][MAX_ATOMS], 
                int   type[MAX_ATOMS], 
                Real vt[MAX_TORS][SPACE],
		        char  rms_ref_crds[MAX_CHARS],
		        Real torsFreeEnergy,
                Boole B_write_all_clusmem,
                int ligand_is_inhibitor,
                Boole B_template,
                Real template_energy[MAX_ATOMS], // template energy value for each atom
                Real template_stddev[MAX_ATOMS], // and standard deviation of this energy
                int   outlev,
                int   ignore_inter[MAX_ATOMS],
                const Boole   B_include_1_4_interactions,
                const Real scale_1_4,

                const ParameterEntry parameterArray[MAX_MAPS],
                const Real unbound_internal_FE,
                
                GridMapSetInfo *info);
#endif
