#ifndef TRILINTERP
#define TRILINTERP

#define ENERGYPENALTY 500.0	/* Energy factor which is multiplied by distance
				   from centre of grid, to penalize atoms
				   outside grid */

#define ALL_ATOMS_INSIDE_GRID 0
#define SOME_ATOMS_OUTSIDE_GRID 1

#define NULL_EVDW ((Real *)NULL)
#define NULL_ELEC ((Real *)NULL)
#define NULL_EVDW_TOTAL ((Real *)NULL)
#define NULL_ELEC_TOTAL ((Real *)NULL)
#define NULL_IGNORE_INTERMOL ((int *)NULL)

#include "constants.h"
#include "structs.h"

Real  trilinterp( CONST_INT first_atom,
 CONST_INT last_atom,
 CONST_FLOAT tcoord[MAX_ATOMS][SPACE], // temporary coordinates
 CONST_FLOAT charge[MAX_ATOMS], // partial atomic charges
 CONST_FLOAT abs_charge[MAX_ATOMS], 
 CONST_INT   type[MAX_ATOMS], // atom type of each atom
 CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],    //  intermolecular interaction energies
 GridMapSetInfo *info, // info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
 int some_atoms_outside_grid, // boolean
 int ignore_inter[MAX_ATOMS], // array of booleans, says to ignore computation intermolecular energies per atom
 Real elec[MAX_ATOMS], // set if not NULL - electrostatic energies, atom by atom
 Real emap[MAX_ATOMS],  // set if not NULL - intermolecular energies
 Real *p_elec_total, // set if not NULL - total electrostatic energy
 Real *p_emap_total // set if not NULL - total intermolecular energy
 );

#endif

