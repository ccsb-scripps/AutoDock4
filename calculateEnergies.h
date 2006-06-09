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
    int ** const         nonbondlist,               // input  list of nonbonds
    const EnergyTables   *ptr_ad_energy_tables,     // input  pointer to AutoDock intermolecular, dielectric, solvation lookup tables
    const int            Nnb,                       // input  total number of nonbonds
    const Boole          B_calcIntElec,             // input  boolean whether we must calculate internal electrostatics
    const Real           q1q2[MAX_NONBONDS],        // input  product of partial charges for each nonbonded pair of atoms
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
