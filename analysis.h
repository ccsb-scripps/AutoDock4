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
                FloatOrDouble charge[MAX_ATOMS], 
	              Boole B_calcIntElec,
	              FloatOrDouble q1q2[MAX_NONBONDS],
                FloatOrDouble clus_rms_tol, 
                FloatOrDouble crdpdb[MAX_ATOMS][SPACE], 
                FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                FloatOrDouble inv_spacing, 
                FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                FloatOrDouble econf[MAX_RUNS], 
                int   irunmax, 
                FloatOrDouble xlo, 
                FloatOrDouble ylo, 
                FloatOrDouble zlo, 
                int   natom, 
                int   nonbondlist[MAX_NONBONDS][4], 
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
                int   outlev);
#endif
