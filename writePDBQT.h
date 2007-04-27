/*

 $Id: writePDBQT.h,v 1.6 2007/04/27 06:01:52 garrett Exp $

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

#ifndef _WRITEPDBQT
#define _WRITEPDBQT

#include "structs.h"
#include "constants.h"
#include "printEnergies.h"
#include "trilinterp.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

void writePDBQT(int irun,FourByteLong seed[2],
                    char  smFileName[MAX_CHARS],
                    char  dpfFN[MAX_CHARS],
                    Real sml_center[SPACE],
                    State state,
                    int   ntor,
                    Real (*Ptr_eintra),
                    Real (*Ptr_einter),
                    int   natom,
                    char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    Real crd[MAX_ATOMS][SPACE],
                    Real emap[MAX_ATOMS],
                    Real elec[MAX_ATOMS],
                    Real charge[MAX_ATOMS],
                    Real abs_charge[MAX_ATOMS],
                    Real qsp_abs_charge[MAX_ATOMS],
                    int ligand_is_inhibitor,
                    Real torsFreeEnergy,
                    Real vt[MAX_TORS][SPACE],
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    Real crdpdb[MAX_ATOMS][SPACE],
                    NonbondParam *nonbondlist,
                    EnergyTables *ptr_ad_energy_tables,
                    int   type[MAX_ATOMS],
                    int   Nnb,
                    Boole B_calcIntElec,
                    Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    int outlev,
                    int   ignore_inter[MAX_ATOMS],
                    const Boole         B_include_1_4_interactions,
                    const Real scale_1_4,
                    const ParameterEntry parameterArray[MAX_MAPS],
                    const Real unbound_internal_FE,
                    GridMapSetInfo *info,
                    int state_type,  // 0 means unbound, 1 means docked
                    char PDBQT_record[MAX_RECORDS][LINE_LEN],
                    Boole B_use_non_bond_cutoff,
                    Boole B_have_flexible_residues
                    );

void print_PDBQT( FILE *logFile, 
                  const int true_ligand_atoms,
                  const char atomstuff[MAX_ATOMS][MAX_CHARS],
                  const Real crdpdb[MAX_ATOMS][SPACE],
                  const Real charge[MAX_ATOMS],
                  const ParameterEntry parameterArray[MAX_MAPS],
                  const int type[MAX_ATOMS],
                  const char prefix[MAX_CHARS] );

#endif
