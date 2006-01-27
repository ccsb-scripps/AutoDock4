#ifndef TRILINTERP
#define TRILINTERP

#define ENERGYPENALTY 500.0	/* Energy factor which is multiplied by distance
				   from centre of grid, to penalize atoms
				   outside grid */

#define ALL_ATOMS_INSIDE_GRID 0
#define SOME_ATOMS_OUTSIDE_GRID 1

#define NULL_EVDW ((FloatOrDouble *)NULL)
#define NULL_ELEC ((FloatOrDouble *)NULL)
#define NULL_EVDW_TOTAL ((FloatOrDouble *)NULL)
#define NULL_ELEC_TOTAL ((FloatOrDouble *)NULL)
#define NULL_IGNORE_INTERMOL ((int *)NULL)

#include "constants.h"
#include "structs.h"

FloatOrDouble  trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], // temporary coordinates
 CONST_FLOAT charge[MAX_ATOMS], // partial atomic charges
 CONST_FLOAT abs_charge[MAX_ATOMS], 
 CONST_INT   type[MAX_ATOMS], // atom type of each atom
 CONST_INT   total_atoms, // number of atoms
 CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],    //  intermolecular interaction energies
 GridMapSetInfo *info, // info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
 int some_atoms_outside_grid, // boolean
 int ignore_inter[MAX_ATOMS], // array of booleans, says to ignore computation intermolecular energies per atom
 FloatOrDouble elec[MAX_ATOMS], // set if not NULL - electrostatic energies, atom by atom
 FloatOrDouble emap[MAX_ATOMS],  // set if not NULL - intermolecular energies
 FloatOrDouble *p_elec_total, // set if not NULL - total electrostatic energy
 FloatOrDouble *p_emap_total // set if not NULL - total intermolecular energy
 );

FloatOrDouble  template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], // temporary coordinates
 CONST_FLOAT charge[MAX_ATOMS], // partial atomic charges
 CONST_FLOAT abs_charge[MAX_ATOMS], 
 CONST_INT   type[MAX_ATOMS], // atom type of each atom
 CONST_INT   total_atoms, // number of atoms
 CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],    //  intermolecular interaction energies
 GridMapSetInfo *info, // info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
 int some_atoms_outside_grid, // boolean
 int ignore_inter[MAX_ATOMS], // array of booleans, says to ignore computation intermolecular energies per atom
 CONST_FLOAT template_energy[MAX_ATOMS],
 CONST_FLOAT template_stddev[MAX_ATOMS],
 FloatOrDouble elec[MAX_ATOMS], // set if not NULL - electrostatic energies, atom by atom
 FloatOrDouble emap[MAX_ATOMS],  // set if not NULL - intermolecular energies
 FloatOrDouble *p_elec_total, // set if not NULL - total electrostatic energy
 FloatOrDouble *p_emap_total // set if not NULL - total intermolecular energy
 );
#endif

