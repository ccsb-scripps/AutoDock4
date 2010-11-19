/*

 $Id: cmdmode.h,v 1.13 2007/04/27 06:01:48 garrett Exp $

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

#ifndef CMDMODE
#define CMDMODE

#include "constants.h"
#include "set_cmd_io_std.h"
#include "print_2x.h"
#include "parse_com_line.h"
#include "strindex.h"
#include "print_avsfld.h"
#include "printEnergies.h"
#include "success.h"
#include "readPDBQT.h"
#include "get_atom_type.h"
#include "timesys.h"
#include "eintcalPrint.h"
#include "trilinterp.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "parse_trj_line.h"
#include "input_state.h"
#include "openfile.h"

int   cmdmode( int natom,
             Clock jobStart,
             struct tms tms_jobStart,
             Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

             EnergyTables *ptr_ad_energy_tables,

             Real WallEnergy,
             Real vt[MAX_TORS][SPACE],
             int   tlist[MAX_TORS][MAX_ATOMS],
             int   ntor,
             int   Nnb,
             NonbondParam *nonbondlist,
             char  atomstuff[MAX_ATOMS][MAX_CHARS],
             Real crdpdb[MAX_ATOMS][SPACE],
             char  hostnm[MAX_CHARS],
             int   type[MAX_ATOMS],
             Real charge[MAX_ATOMS],
             Real abs_charge[MAX_ATOMS],
             Real qsp_abs_charge[MAX_ATOMS],
             Boole B_calcIntElec,
             char  atm_typ_str[ATOM_MAPS],
             Real torsFreeEnergy,
             int ligand_is_inhibitor,
             int ignore_inter[MAX_ATOMS],
             const Boole         B_include_1_4_interactions,
             const Real scale_1_4,

             const ParameterEntry parameterArray[MAX_MAPS],
             const Real unbound_internal_FE,

             GridMapSetInfo *info,
             Boole B_have_flexible_residues,
             Boole B_use_non_bond_cutoff
             );
#endif
