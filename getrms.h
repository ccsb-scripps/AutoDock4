#ifndef GETRMS
#define GETRMS

#include "constants.h"

FloatOrDouble  getrms( FloatOrDouble Crd[MAX_ATOMS][SPACE], 
               FloatOrDouble CrdRef[MAX_ATOMS][SPACE], 
               Boole B_symmetry_flag, 
               int   natom, 
               int   type[MAX_ATOMS] );
#endif
