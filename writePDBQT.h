/*

 $Id: writePDBQT.h,v 1.13 2010/08/27 00:05:09 mp Exp $

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

#ifndef _WRITEPDBQT
#define _WRITEPDBQT

#include "structs.h"
#include "constants.h"
#include "printEnergies.h"
#include "trilinterp.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

void writePDBQT(const int irun, const FourByteLong seed[2],
                    const char  *const smFileName,
                    const char  *dpfFN,
                    const Real sml_center[SPACE],
                    /* not const */ State state,
                    const int   ntor,
                    /* not const */ Real (*const Ptr_eintra),
                    /* not const */ Real (*const Ptr_einter),
                    const int   natom,
                    const char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    /* not const */ Real crd[MAX_ATOMS][SPACE],
                    /* not const */ Real emap[MAX_ATOMS],
                    /* not const */ Real elec[MAX_ATOMS],
                    const Real charge[MAX_ATOMS],
                    const Real abs_charge[MAX_ATOMS],
                    const Real qsp_abs_charge[MAX_ATOMS],
                    const int ligand_is_inhibitor,
                    const Real torsFreeEnergy,
                    const Real vt[MAX_TORS][SPACE],
                    const int   tlist[MAX_TORS][MAX_ATOMS],
                    const Real crdpdb[MAX_ATOMS][SPACE],
                    const NonbondParam *const nonbondlist,
                    const EnergyTables *const ptr_ad_energy_tables,
                    const int   type[MAX_ATOMS],
                    const int   Nnb,
                    const Boole B_calcIntElec,
                #include "map_declare.h"
                    const int outlev,
                    const int   ignore_inter[MAX_ATOMS],
                    const Boole B_include_1_4_interactions,
                    const Real scale_1_4,
                    const ParameterEntry parameterArray[MAX_ATOM_TYPES], // input  nonbond and desolvation parameters
                    const Real unbound_internal_FE,
                    const GridMapSetInfo *const info,
                    const int state_type,  // 0 means unbound, 1 means docked
                    const char PDBQT_record[MAX_RECORDS][LINE_LEN],
                    const Boole B_use_non_bond_cutoff,
                    const Boole B_have_flexible_residues,
                    const Unbound_Model ad4_unbound_model
                    );

void print_PDBQT( FILE *const logFile, 
                  const int true_ligand_atoms,
                  const char atomstuff[MAX_ATOMS][MAX_CHARS],
                  const Real crdpdb[MAX_ATOMS][SPACE],
                  const Real charge[MAX_ATOMS],
                  const ParameterEntry parameterArray[MAX_ATOM_TYPES], // input  nonbond and desolvation parameters
                  const int type[MAX_ATOMS],
                  const char prefix[MAX_CHARS] );



void print_PDBQ_atom_resstr( FILE *const logFile, 
                  const char prefix[MAX_CHARS],
                  int atom_num, // 0-origin 
                  const char atomstuff[],
                  const Real crdpdb[MAX_ATOMS][SPACE],
                  const Real vdW,
                  const Real Elec,
                  const Real charge,
                  const char *const suffix //newline or empty
                  );
        

void print_PDBQ_atom_resnum( FILE *const logFile, 
                  const char prefix[MAX_CHARS],
                  const int atom_num, // 0-origin 
                  const char atomstuff[],
                  const int resnum,
                  const Real crdpdb[MAX_ATOMS][SPACE],
                  const Real vdW,
                  const Real Elec,
                  const Real charge,
                  const char *const suffix //newline or empty
                  );
        
#endif
