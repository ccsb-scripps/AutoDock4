#ifndef EVALUATE_ENERGY
#define EVALUATE_ENERGY

#include <stdio.h>
#include "constants.h"
#include "trilinterp.h"
#include "eintcal.h"

FloatOrDouble evaluate_energy( 
    FloatOrDouble crd[MAX_ATOMS][SPACE],
    FloatOrDouble charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    FloatOrDouble inv_spacing,
    FloatOrDouble xlo,
    FloatOrDouble ylo,
    FloatOrDouble zlo,
    int   nonbondlist[MAX_NONBONDS][2],
    FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
    int   Nnb,
    Boole B_calcIntElec,
    FloatOrDouble q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    Boole B_isTorConstrained[MAX_TORS],
    State now,
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS]);

/*
FloatOrDouble evaluate_energy( 
    FloatOrDouble crd[MAX_ATOMS][SPACE],
    FloatOrDouble charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    FloatOrDouble inv_spacing,
    FloatOrDouble xlo,
    FloatOrDouble ylo,
    FloatOrDouble zlo,
/ *    FloatOrDouble *Addr_eintra, * /
    int   nonbondlist[MAX_NONBONDS][2],
    FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
    int   Nnb,
    Boole B_calcIntElec,
    FloatOrDouble q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    / * int   ntor, * /
    Boole B_isTorConstrained[MAX_TORS],
    State now, / * WAS: FloatOrDouble torNow[MAX_TORS], * /
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS]);
    */

#endif
