/*

 $Id: simanneal.h,v 1.19 2010/08/27 00:05:08 mp Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

void simanneal( int   *const P_nconf, 
                const int   Nnb, 
                const Real WallEnergy, 
                const char  atomstuff[MAX_ATOMS][MAX_CHARS], 
                Real charge[MAX_ATOMS], 
                Real abs_charge[MAX_ATOMS], 
                Real qsp_abs_charge[MAX_ATOMS], 
                const Boole B_calcIntElec,
                Real crd[MAX_ATOMS][SPACE], 
                Real crdpdb[MAX_ATOMS][SPACE], 
                const char  *const dpfFN,
                
                    EnergyTables *const ptr_ad_energy_tables,

                /* not const */ Real econf[MAX_RUNS], 
                const Boole B_either, 
                Real elec[MAX_ATOMS], 
                Real emap[MAX_ATOMS], 
                const int   icyclemax, 
                const int   irunmax, 
                const Clock jobStart, 
                #include "map_declare.h"
                const int   naccmax, 
                const int   natom, 
                NonbondParam *const nonbondlist, 
                const int   nrejmax, 
                const int   ntor1, 
                const int   ntor, 
                const int   outlev, 

                State sInit, 
		/* not const */ State sHist[MAX_RUNS],        
                const Real qtwFac, 
                const Boole B_qtwReduc, 
                const Real qtwStep0, 
                const Boole B_selectmin, 
                const char  *const smFileName,
                const Real sml_center[SPACE], 
                const Real RT0, 
                const Boole B_RTChange, 
                const Real RTFac, 
                struct tms tms_jobStart, 
                const int   tlist[MAX_TORS][MAX_ATOMS], 
                const Real torFac, 
                const Boole B_torReduc, 
                const Real torStep0, 
                const char  *const trjFileName,
                const int   trj_cyc_max, 
                const int   trj_cyc_min, 
                const int   trj_freq, 
                const Real trnFac, 
                const Boole B_trnReduc, 
                const Real trnStep0, 
                const int   type[MAX_ATOMS], 
                Real vt[MAX_TORS][SPACE], 
                const Boole B_write_trj, 
                const Boole B_constrain, 
                const int   atomC1, 
                const int   atomC2, 
                const Real sqlower, 
                const Real squpper,
                const Boole B_linear_schedule,
                const Real RTreduc,
                /*Real maxrad,*/
                const Boole B_watch,
                const char  *const FN_watch,
                const Boole B_isGaussTorCon,
                /* not const */ unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                const Boole B_isTorConstrained[MAX_TORS],
                const Boole B_ShowTorE,
                /* not const */ unsigned short US_TorE[MAX_TORS],
                const Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                const int   N_con[MAX_TORS],
                const Boole B_RandomTran0,
                const Boole B_RandomQuat0,
                const Boole B_RandomDihe0,
                const Real e0max,
                const Real torsFreeEnergy,
                const int   MaxRetries,
                const int   ligand_is_inhibitor,
                const int   ignore_inter[MAX_ATOMS],
                const Boole         B_include_1_4_interactions,
                const Real scale_1_4,
                const Real scale_eintermol,
                const ParameterEntry parameterArray[MAX_ATOM_TYPES], // input  nonbond and desolvation parameters
                const Real unbound_internal_FE,

                GridMapSetInfo *const info,
                const Boole B_use_non_bond_cutoff,
                const Boole B_have_flexible_residues,
                const char  PDBQT_record[MAX_RECORDS][LINE_LEN],
                const Unbound_Model ad4_unbound_model);

#endif
