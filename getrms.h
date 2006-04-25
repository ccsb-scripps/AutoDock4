#ifndef GETRMS
#define GETRMS

#include "constants.h"

Real  getrms( Real Crd[MAX_ATOMS][SPACE], 
               Real CrdRef[MAX_ATOMS][SPACE], 
               Boole B_symmetry_flag, 
               int   natom, 
               int   type[MAX_ATOMS] );
#endif
