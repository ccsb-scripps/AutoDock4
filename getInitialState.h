/*

 $Id: getInitialState.h,v 1.12 2007/04/27 06:01:48 garrett Exp $

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

#ifndef GETINITIALSTATE
#define GETINITIALSTATE

#include "constants.h"
#include "qmultiply.h"
#include "stateLibrary.h"
#include "initautodock.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "cnv_state_to_coords.h"
#include "prInitialState.h"
#include "timesys.h"

void getInitialState(  
            Real *Addr_e0,
            Real e0max,

	    State *sInit,
	    State *sMin,
	    State *sLast,

            Boole B_RandomTran0,
            Boole B_RandomQuat0,
            Boole B_RandomDihe0,

            Real charge[MAX_ATOMS],
            Real abs_charge[MAX_ATOMS],
            Real qsp_abs_charge[MAX_ATOMS],
            Real crd[MAX_ATOMS][SPACE],
            Real crdpdb[MAX_ATOMS][SPACE],
            char  atomstuff[MAX_ATOMS][MAX_CHARS],
            Real elec[MAX_ATOMS],
            Real emap[MAX_ATOMS],

            EnergyTables *ptr_ad_energy_tables,

            Boole B_calcIntElec,
            Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int   natom,
            int   Nnb,
            NonbondParam *nonbondlist,
            int   ntor,
            int   tlist[MAX_TORS][MAX_ATOMS],
            int   type[MAX_ATOMS],
            Real vt[MAX_TORS][SPACE],
            int   irun1,
            int   outlev,
	        int   MaxRetries,

	        Real torsFreeEnergy,

            int   ligand_is_inhibitor,

            int ignore_inter[MAX_ATOMS],

            const Boole         B_include_1_4_interactions,
            const Real scale_1_4,

            const ParameterEntry parameterArray[MAX_MAPS],

            const Real unbound_internal_FE,

            GridMapSetInfo *info,
            Boole B_use_non_bond_cutoff,
            Boole B_have_flexible_residues

           );

#endif
