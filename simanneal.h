
#ifndef SIMANNEAL
#define SIMANNEAL
#include "constants.h"
#include "print_2x.h"
#include "timesys.h"
#include "getInitialState.h"
#include "mkNewState.h"
#include "stateLibrary.h"
#include "output_state.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "writePDBQ.h"
void simanneal( int   *P_nconf, 
                int   Nnb, 
                FloatOrDouble WallEnergy, 
                char  atomstuff[MAX_ATOMS][MAX_CHARS], 
                FloatOrDouble charge[MAX_ATOMS], 
		Boole B_calcIntElec,
		FloatOrDouble q1q2[MAX_NONBONDS],
                FloatOrDouble crd[MAX_ATOMS][SPACE], 
                FloatOrDouble crdpdb[MAX_ATOMS][SPACE], 
		char  dpfFN[MAX_CHARS],
                FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                FloatOrDouble econf[MAX_RUNS], 
                Boole B_either, 
                FloatOrDouble elec[MAX_ATOMS], 
                FloatOrDouble emap[MAX_ATOMS], 
                FloatOrDouble xhi, 
                FloatOrDouble yhi, 
                FloatOrDouble zhi, 
                int   icyclemax, 
                FloatOrDouble inv_spacing, 
                int   irunmax, 
                Clock jobStart, 
                FloatOrDouble xlo, 
                FloatOrDouble ylo, 
                FloatOrDouble zlo, 
                FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                int   naccmax, 
                int   natom, 
                int   nonbondlist[MAX_NONBONDS][4], 
                int   nrejmax, 
                int   ntor1, 
                int   ntor, 
                int   outlev, 

                State sInit,  		/* FloatOrDouble qtn0[QUAT], 
				           FloatOrDouble tor0[MAX_TORS], */
                State sHist[MAX_RUNS],	/* FloatOrDouble qtnHist[MAX_RUNS][QUAT], 
				 	   FloatOrDouble torHist[MAX_RUNS][MAX_TORS], */
                FloatOrDouble qtwFac, 
                Boole B_qtwReduc, 
                FloatOrDouble qtwStep0, 
                Boole B_selectmin, 
                char  smFileName[MAX_CHARS], 
                FloatOrDouble sml_center[SPACE], 
                FloatOrDouble RT0, 
                Boole B_RTChange, 
                FloatOrDouble RTFac, 
                struct tms tms_jobStart, 
                int   tlist[MAX_TORS][MAX_ATOMS], 
                FloatOrDouble torFac, 
                Boole B_torReduc, 
                FloatOrDouble torStep0, 
                char  trjFileName[MAX_CHARS], 
                int   trj_cyc_max, 
                int   trj_cyc_min, 
                int   trj_freq, 
                FloatOrDouble trnFac, 
                Boole B_trnReduc, 
                FloatOrDouble trnStep0, 
                int   type[MAX_ATOMS], 
                FloatOrDouble vt[MAX_TORS][SPACE], 
                Boole B_write_trj, 
                Boole B_constrain, 
                int   atomC1, 
                int   atomC2, 
                FloatOrDouble sqlower, 
                FloatOrDouble squpper,
		Boole B_linear_schedule,
		FloatOrDouble RTreduc,
		/*FloatOrDouble maxrad,*/
		Boole B_watch,
		char  FN_watch[MAX_CHARS],
		Boole B_isGaussTorCon,
		unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		Boole B_isTorConstrained[MAX_TORS],
		Boole B_ShowTorE,
		unsigned short US_TorE[MAX_TORS],
	        FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		int   N_con[MAX_TORS],
                Boole B_RandomTran0,
                Boole B_RandomQuat0,
                Boole B_RandomDihe0,
                FloatOrDouble e0max,
                FloatOrDouble torsFreeEnergy,
                int   MaxRetries,
                int   ligand_is_inhibitor,
                int   ignore_inter[MAX_ATOMS]);

#endif
