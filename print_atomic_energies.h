
#ifndef PRINT_ATOMIC_ENERGIES
#define PRINT_ATOMIC_ENERGIES
#include "constants.h"
void  print_atomic_energies(int   natom,
                            char  atomstuff[MAX_ATOMS][MAX_CHARS],
                            int   type[MAX_ATOMS],
                            FloatOrDouble emap[MAX_ATOMS],
                            FloatOrDouble elec[MAX_ATOMS],
                            FloatOrDouble charge[MAX_ATOMS] );
#endif
