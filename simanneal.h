
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
                float WallEnergy, 
                char  atomstuff[MAX_ATOMS][MAX_CHARS], 
                float charge[MAX_ATOMS], 
		Boole B_calcIntElec,
		float q1q2[MAX_NONBONDS],
                float crd[MAX_ATOMS][SPACE], 
                float crdpdb[MAX_ATOMS][SPACE], 
		char  dpfFN[MAX_CHARS],
                float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                float econf[MAX_RUNS], 
                Boole B_either, 
                float elec[MAX_ATOMS], 
                float emap[MAX_ATOMS], 
                float xhi, 
                float yhi, 
                float zhi, 
                int   icyclemax, 
                float inv_spacing, 
                int   irunmax, 
                Clock jobStart, 
                float xlo, 
                float ylo, 
                float zlo, 
                float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                int   naccmax, 
                int   natom, 
                int   nonbondlist[MAX_NONBONDS][2], 
                int   nrejmax, 
                int   ntor1, 
                int   ntor, 
                int   outlev, 

                State sInit,  		/* float qtn0[QUAT], 
				           float tor0[MAX_TORS], */
                State sHist[MAX_RUNS],	/* float qtnHist[MAX_RUNS][QUAT], 
				 	   float torHist[MAX_RUNS][MAX_TORS], */
                float qtwFac, 
                Boole B_qtwReduc, 
                float qtwStep0, 
                Boole B_selectmin, 
                char  smFileName[MAX_CHARS], 
                float sml_center[SPACE], 
                float RT0, 
                Boole B_RTChange, 
                float RTFac, 
                struct tms tms_jobStart, 
                int   tlist[MAX_TORS][MAX_ATOMS], 
                float torFac, 
                Boole B_torReduc, 
                float torStep0, 
                char  trjFileName[MAX_CHARS], 
                int   trj_cyc_max, 
                int   trj_cyc_min, 
                int   trj_freq, 
                float trnFac, 
                Boole B_trnReduc, 
                float trnStep0, 
                int   type[MAX_ATOMS], 
                float vt[MAX_TORS][SPACE], 
                Boole B_write_trj, 
                Boole B_constrain, 
                int   atomC1, 
                int   atomC2, 
                float sqlower, 
                float squpper,
		Boole B_linear_schedule,
		float RTreduc,
		/*float maxrad,*/
		Boole B_watch,
		char  FN_watch[MAX_CHARS],
		Boole B_isGaussTorCon,
		unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		Boole B_isTorConstrained[MAX_TORS],
		Boole B_ShowTorE,
		unsigned short US_TorE[MAX_TORS],
	        float F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		int   N_con[MAX_TORS],
                Boole B_RandomTran0,
                Boole B_RandomQuat0,
                Boole B_RandomDihe0,
                float e0max,
                float torsFreeEnergy,
                int   MaxRetries,
                int   ligand_is_inhibitor);

#endif
