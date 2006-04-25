
#ifndef PRINT_ATOMIC_ENERGIES
#define PRINT_ATOMIC_ENERGIES
#include "constants.h"
void  print_atomic_energies(int   natom,
                            char  atomstuff[MAX_ATOMS][MAX_CHARS],
                            int   type[MAX_ATOMS],
                            Real emap[MAX_ATOMS],
                            Real elec[MAX_ATOMS],
                            Real charge[MAX_ATOMS] );
#endif
