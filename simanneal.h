/*

 $Id: simanneal.h,v 1.12 2007/04/27 06:01:51 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

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
#include "writePDBQT.h"

void simanneal( int   *P_nconf, 
                int   Nnb, 
                Real WallEnergy, 
                char  atomstuff[MAX_ATOMS][MAX_CHARS], 
                Real charge[MAX_ATOMS], 
                Real abs_charge[MAX_ATOMS], 
                Real qsp_abs_charge[MAX_ATOMS], 
                Boole B_calcIntElec,
                Real crd[MAX_ATOMS][SPACE], 
                Real crdpdb[MAX_ATOMS][SPACE], 
                char  dpfFN[MAX_CHARS],
                
                    EnergyTables *ptr_ad_energy_tables,

                Real econf[MAX_RUNS], 
                Boole B_either, 
                Real elec[MAX_ATOMS], 
                Real emap[MAX_ATOMS], 
                int   icyclemax, 
                int   irunmax, 
                Clock jobStart, 
                Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                int   naccmax, 
                int   natom, 
                NonbondParam *nonbondlist, 
                int   nrejmax, 
                int   ntor1, 
                int   ntor, 
                int   outlev, 

                State sInit,                  /* Real qtn0[QUAT], 
                                           Real tor0[MAX_TORS], */
                State sHist[MAX_RUNS],        /* Real qtnHist[MAX_RUNS][QUAT], 
                                            Real torHist[MAX_RUNS][MAX_TORS], */
                Real qtwFac, 
                Boole B_qtwReduc, 
                Real qtwStep0, 
                Boole B_selectmin, 
                char  smFileName[MAX_CHARS], 
                Real sml_center[SPACE], 
                Real RT0, 
                Boole B_RTChange, 
                Real RTFac, 
                struct tms tms_jobStart, 
                int   tlist[MAX_TORS][MAX_ATOMS], 
                Real torFac, 
                Boole B_torReduc, 
                Real torStep0, 
                char  trjFileName[MAX_CHARS], 
                int   trj_cyc_max, 
                int   trj_cyc_min, 
                int   trj_freq, 
                Real trnFac, 
                Boole B_trnReduc, 
                Real trnStep0, 
                int   type[MAX_ATOMS], 
                Real vt[MAX_TORS][SPACE], 
                Boole B_write_trj, 
                Boole B_constrain, 
                int   atomC1, 
                int   atomC2, 
                Real sqlower, 
                Real squpper,
                Boole B_linear_schedule,
                Real RTreduc,
                /*Real maxrad,*/
                Boole B_watch,
                char  FN_watch[MAX_CHARS],
                Boole B_isGaussTorCon,
                unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                Boole B_isTorConstrained[MAX_TORS],
                Boole B_ShowTorE,
                unsigned short US_TorE[MAX_TORS],
                Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                int   N_con[MAX_TORS],
                Boole B_RandomTran0,
                Boole B_RandomQuat0,
                Boole B_RandomDihe0,
                Real e0max,
                Real torsFreeEnergy,
                int   MaxRetries,
                int   ligand_is_inhibitor,
                int   ignore_inter[MAX_ATOMS],
                const Boole         B_include_1_4_interactions,
                const Real scale_1_4,
                const ParameterEntry parameterArray[MAX_MAPS],
                const Real unbound_internal_FE,

                GridMapSetInfo *info,
                Boole B_use_non_bond_cutoff,
                Boole B_have_flexible_residues,
                char  PDBQT_record[MAX_RECORDS][LINE_LEN]);

#endif
