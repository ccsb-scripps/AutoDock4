/*

 $Id: calculateEnergies.cc,v 1.3 2007/04/27 06:01:47 garrett Exp $

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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* calculateEnergies.cc */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "calculateEnergies.h"
#include "trilinterp.h"
#include "eintcal.h"

extern FILE *logFile;
extern int true_ligand_atoms;
extern int Nnb_array[3];
extern Real nb_group_energy[3];

// EnergyBreakdown eb;
// eb = calculateEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues,
//      tcoord, charge, abs_charge, type, map, ptr_info, B_outside, 
//      ignore_inter, elec, emap, p_elec_total, p_emap_total,
//      nonbondlist, ptr_ad_energy_tables, Nnb, B_calcIntElec, 
//      B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff );

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
    const Boole          B_use_non_bond_cutoff     // input  boolean whether to use a nonbond distance cutoff

)

{
    EnergyBreakdown eb;

    initialise_energy_breakdown( &eb, torsFreeEnergy, unbound_internal_FE );

    // computing trilinear-interpolated energies from atom = 0 to atom < true_ligand_atoms
    // gives the intermolecular energy between the ligand and the protein
    eb.e_inter_moving_fixed = trilinterp( 0, true_ligand_atoms, tcoord, charge, abs_charge, type, map, 
                         info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                         ignore_inter, elec, emap,
                         p_elec_total, p_emap_total);

    if (B_have_flexible_residues) {
        // computing trilinear-interpolated energies from atom = true_ligand_atoms to atom < true_ligand_atoms
        // gives the intramolecular energy within the protein
        // we can ignore the elec_total and emap_total breakdown here
        eb.e_intra_moving_fixed_rec = trilinterp( true_ligand_atoms, natom, tcoord, charge, abs_charge, type, map, 
                             info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                             ignore_inter, elec, emap,
                             NULL, NULL);
    }

    if (ntor > 0) {
        // computing all the nonbond interaction energies fills nb_group_energy[3] array
        // with intramolecular energy of ligand, intermolecular energy, and intramolecular energy of receptor
        (void) eintcal( nonbondlist, ptr_ad_energy_tables, tcoord, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues ) ;
        
        eb.e_intra_moving_moving_lig = nb_group_energy[INTRA_LIGAND];
        eb.e_inter_moving_moving = nb_group_energy[INTER];
        eb.e_intra_moving_moving_rec = nb_group_energy[INTRA_RECEPTOR];
    }

    // update the totals in the energy breakdown structure
    update_energy_breakdown( &eb );

    return eb;
} // calculateEnergies()

void update_energy_breakdown( EnergyBreakdown * eb )
{
    // total intermolecular energy = (1) + (4)
    eb->e_inter     = eb->e_inter_moving_fixed + eb->e_inter_moving_moving;

    // total intramolecular energy = (3)  +  (2) + (5)
    eb->e_intra     = eb->e_intra_moving_moving_lig + eb->e_intra_moving_fixed_rec + eb->e_intra_moving_moving_rec;

    // ligand intramolecular energy = (3)
    eb->e_intra_lig = eb->e_intra_moving_moving_lig;

    // receptor intramolecular energy = (2) + (5)
    eb->e_intra_rec = eb->e_intra_moving_fixed_rec + eb->e_intra_moving_moving_rec;

    // estimated free energy upon binding
    eb->deltaG = eb->e_inter + eb->e_intra + eb->e_torsFreeEnergy - eb->e_unbound_internal_FE;
}

void initialise_energy_breakdown ( EnergyBreakdown * eb,
                                   Real torsFreeEnergy, 
                                   Real unbound_internal_FE )
{
    eb->e_inter_moving_fixed = 0.0;      // (1)  // trilinterp( 0, true_ligand_atoms, ...)
    eb->e_intra_moving_fixed_rec = 0.0;  // (2)  // trilinterp( true_ligand_atoms, natom, ...)
    eb->e_intra_moving_moving_lig = 0.0; // (3)  // eintcal( 0, nb_array[0], ...)            // nb_group_energy[INTRA_LIGAND]
    eb->e_inter_moving_moving = 0.0;     // (4)  // eintcal( nb_array[0], nb_array[1], ...)  // nb_group_energy[INTER]
    eb->e_intra_moving_moving_rec = 0.0; // (5)  // eintcal( nb_array[1], nb_array[2], ...)  // nb_group_energy[INTRA_RECEPTOR]

    eb->e_inter = 0.0;                   // total    intermolecular energy = (1) + (4)
    eb->e_intra = 0.0;                   // total    intramolecular energy = (3)  +  (2) + (5)
    eb->e_intra_lig = 0.0;               // ligand   intramolecular energy = (3)
    eb->e_intra_rec = 0.0;               // receptor intramolecular energy = (2) + (5)

    eb->e_torsFreeEnergy = torsFreeEnergy;            // empirical torsional free energy penalty
    eb->e_unbound_internal_FE = unbound_internal_FE;  // computed internal free energy of the unbound state
    eb->deltaG = 0.0;                    // estimated change in free energy upon binding
}

// EOF
