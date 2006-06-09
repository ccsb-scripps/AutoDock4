#ifndef PRINITIALSTATE
#define PRINITIALSTATE

#include "constants.h"
#include "print_atomic_energies.h"
#include "printEnergies.h"

void  prInitialState( 
    EnergyBreakdown *eb,
    int natom,
    Real crd[MAX_ATOMS][SPACE],
    char atomstuff[MAX_ATOMS][MAX_CHARS],
    int type[MAX_ATOMS],
    Real emap[MAX_ATOMS],
    Real elec[MAX_ATOMS],
    Real charge[MAX_ATOMS],
    int ligand_is_inhibitor,
    Boole B_have_flexible_residues );

#endif
