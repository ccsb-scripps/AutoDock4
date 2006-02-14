/* investigate.h */


#ifndef INVESTIGATE
#define INVESTIGATE

#include "constants.h"
#include "getpdbcrds.h"
#include "mkRandomState.h"
#include "cnv_state_to_coords.h"
#include "getrms.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "changeState.h"
#include "stateLibrary.h"

#define NUMRMSBINS 80 /* int   NumRmsBins = 40; // NumRmsBins = MaxRms / RmsBinSize; */

extern FILE *logFile;
extern char *programname;


void investigate(
                int   Nnb,
                FloatOrDouble charge[MAX_ATOMS],
                FloatOrDouble abs_charge[MAX_ATOMS],
                FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                Boole B_calcIntElec,
                FloatOrDouble q1q2[MAX_NONBONDS],
                FloatOrDouble crd[MAX_ATOMS][SPACE],
                FloatOrDouble crdpdb[MAX_ATOMS][SPACE],

                EnergyTables *ptr_ad_energy_tables,

                int   maxTests,
                FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                int   natom,
                int   **nonbondlist,
                int   ntor,
                int   outlev,
                int   tlist[MAX_TORS][MAX_ATOMS],
                int   type[MAX_ATOMS],
                FloatOrDouble vt[MAX_TORS][SPACE],
                Boole B_isGaussTorCon,
               unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                Boole B_isTorConstrained[MAX_TORS],
                Boole B_ShowTorE,
               unsigned short US_TorE[MAX_TORS],
                FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                int   N_con[MAX_TORS],
                Boole B_symmetry_flag,
                char  FN_rms_ref_crds[MAX_CHARS],
                int   OutputEveryNTests,
                int   NumLocalTests,
                FloatOrDouble trnStep,
                FloatOrDouble torStep,
                
                int   ignore_inter[MAX_ATOMS],
                
                const Boole         B_include_1_4_interactions,
                const FloatOrDouble scale_1_4,
                
                const ParameterEntry parameterArray[MAX_MAPS],

                const FloatOrDouble unbound_internal_FE,
                GridMapSetInfo *info );
#endif
