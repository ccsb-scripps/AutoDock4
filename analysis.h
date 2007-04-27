/*

 $Id: analysis.h,v 1.15 2007/04/27 06:01:47 garrett Exp $

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
                Real clus_rms_tol, 
                Real crdpdb[MAX_ATOMS][SPACE], 
                
                const EnergyTables *ptr_ad_energy_tables,

                Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                Real econf[MAX_RUNS], 
                int   irunmax, 
                int   natom, 
                NonbondParam *nonbondlist, 
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
                int   outlev,
                int   ignore_inter[MAX_ATOMS],
                const Boole   B_include_1_4_interactions,
                const Real scale_1_4,

                const ParameterEntry parameterArray[MAX_MAPS],
                const Real unbound_internal_FE,
                
                GridMapSetInfo *info,
                Boole B_use_non_bond_cutoff,
                Boole B_have_flexible_residues,
                Boole B_rms_atoms_ligand_only
               );
#endif
