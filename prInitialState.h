
#ifndef PRINITIALSTATE
#define PRINITIALSTATE
#include "constants.h"
#include "print_atomic_energies.h"
#include "printEnergies.h"
void  prInitialState( 
    float einter,
    float eintra,
    float torsFreeEnergy,
    int natom,
    float crd[MAX_ATOMS][SPACE],
    char  atomstuff[MAX_ATOMS][MAX_CHARS],
    int type[MAX_ATOMS],
    float emap[MAX_ATOMS],
    float elec[MAX_ATOMS],
    float charge[MAX_ATOMS],
    int ligand_is_inhibitor);
#endif
