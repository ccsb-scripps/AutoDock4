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
                Real charge[MAX_ATOMS],
                Real abs_charge[MAX_ATOMS],
                Real qsp_abs_charge[MAX_ATOMS],
                Boole B_calcIntElec,
                Real crd[MAX_ATOMS][SPACE],
                Real crdpdb[MAX_ATOMS][SPACE],

                EnergyTables *ptr_ad_energy_tables,

                int   maxTests,
                Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                int   natom,
                NonbondParam *nonbondlist,
                int   ntor,
                int   outlev,
                int   tlist[MAX_TORS][MAX_ATOMS],
                int   type[MAX_ATOMS],
                Real vt[MAX_TORS][SPACE],
                Boole B_isGaussTorCon,
               unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                Boole B_isTorConstrained[MAX_TORS],
                Boole B_ShowTorE,
               unsigned short US_TorE[MAX_TORS],
                Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                int   N_con[MAX_TORS],
                Boole B_symmetry_flag,
                char  FN_rms_ref_crds[MAX_CHARS],
                int   OutputEveryNTests,
                int   NumLocalTests,
                Real trnStep,
                Real torStep,
                
                int   ignore_inter[MAX_ATOMS],
                
                const Boole         B_include_1_4_interactions,
                const Real scale_1_4,
                
                const ParameterEntry parameterArray[MAX_MAPS],

                const Real unbound_internal_FE,
                GridMapSetInfo *info,
                Boole B_use_non_bond_cutoff,
                Boole B_have_flexible_residues);
#endif
