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
		Boole B_calcIntElec,
		FloatOrDouble q1q2[MAX_NONBONDS],
		FloatOrDouble crd[MAX_ATOMS][SPACE],
		FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
		FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
		FloatOrDouble xhi,
		FloatOrDouble yhi,
		FloatOrDouble zhi,
		FloatOrDouble inv_spacing,
		int   maxTests,
		FloatOrDouble xlo,
		FloatOrDouble ylo,
		FloatOrDouble zlo,
		FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
		int   natom,
		int   nonbondlist[MAX_NONBONDS][4],
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
                FloatOrDouble torStep);
#endif
