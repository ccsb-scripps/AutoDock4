
#ifndef PRINITIALSTATE
#define PRINITIALSTATE
#include "constants.h"
#include "print_atomic_energies.h"
#include "printEnergies.h"
void  prInitialState( 
    FloatOrDouble einter,
    FloatOrDouble eintra,
    FloatOrDouble torsFreeEnergy,
    int natom,
    FloatOrDouble crd[MAX_ATOMS][SPACE],
    char  atomstuff[MAX_ATOMS][MAX_CHARS],
    int type[MAX_ATOMS],
    FloatOrDouble emap[MAX_ATOMS],
    FloatOrDouble elec[MAX_ATOMS],
    FloatOrDouble charge[MAX_ATOMS],
    int ligand_is_inhibitor);
#endif
