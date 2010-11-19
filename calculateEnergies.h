/*

 $Id: calculateEnergies.h,v 1.3 2007/04/27 06:01:47 garrett Exp $

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

#ifndef CALCULATEENERGIES
#define CALCULATEENERGIES
#include <stdio.h>
#include "autocomm.h"
#include "constants.h"
#include "structs.h"

EnergyBreakdown calculateEnergies(
    int                  natom,                     // input  number of atoms
    int                  ntor,                      // input  number of torsions
    Real                 unbound_internal_FE,       // input  pre-calculated internal energy of unbound state
    Real                 torsFreeEnergy,            // input  constant times number of freely-rotatable bonds
    Boole                B_have_flexible_residues,  // input  boolean whether we have flexible residues in protein

    // trilinterp
    const Real           tcoord[MAX_ATOMS][SPACE],  // input  coordinates of atoms to be trilinearly-interpolated
    CONST_FLOAT          charge[MAX_ATOMS],         // input  partial atomic charges
    CONST_FLOAT          abs_charge[MAX_ATOMS],     // input  absolute magnitude of partial charges
    CONST_INT            type[MAX_ATOMS],           // input  atom type of each atom
    CONST_FLOAT          map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],    // input  intermolecular interaction energies
    GridMapSetInfo       *info,                     // input  info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
    int                  B_outside,                 // input  boolean whether some atoms are outside grid box
    int                  ignore_inter[MAX_ATOMS],   // input  array of booleans, says to ignore computation intermolecular energies per atom
    Real                 elec[MAX_ATOMS],           // output if not NULL - electrostatic energies, atom by atom
    Real                 emap[MAX_ATOMS],           // output if not NULL - intermolecular energies
    Real                 *p_elec_total,             // output if not NULL - total electrostatic energy
    Real                 *p_emap_total,             // output if not NULL - total intermolecular energy

    // eintcal
    NonbondParam * const         nonbondlist,       // input  list of nonbonds
    const EnergyTables   *ptr_ad_energy_tables,     // input  pointer to AutoDock intermolecular, dielectric, solvation lookup tables
    const int            Nnb,                       // input  total number of nonbonds
    const Boole          B_calcIntElec,             // input  boolean whether we must calculate internal electrostatics
    const Boole          B_include_1_4_interactions,// input  boolean whether to include 1,4 interactions as non-bonds
    const Real           scale_1_4,                 // input  scaling factor for 1,4 interactions, if included
    const Real           qsp_abs_charge[MAX_ATOMS], // input  q-solvation parameters
    const ParameterEntry parameterArray[MAX_MAPS],  // input  nonbond and desolvation parameters
    const Boole          B_use_non_bond_cutoff      // input  boolean whether to use a nonbond distance cutoff

);

void update_energy_breakdown( EnergyBreakdown * eb );

void initialise_energy_breakdown ( EnergyBreakdown * eb,
                                   Real torsFreeEnergy, 
                                   Real unbound_internal_FE );
#endif
