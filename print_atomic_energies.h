
#ifndef PRINT_ATOMIC_ENERGIES
#define PRINT_ATOMIC_ENERGIES
#include "constants.h"
void  print_atomic_energies(int   natom,
                            char  atomstuff[MAX_ATOMS][MAX_CHARS],
                            int   type[MAX_ATOMS],
                            float emap[MAX_ATOMS],
                            float elec[MAX_ATOMS],
                            float charge[MAX_ATOMS] );
#endif
