#ifndef GETRMS
#define GETRMS

#include "constants.h"

float  getrms( float Crd[MAX_ATOMS][SPACE], 
               float CrdRef[MAX_ATOMS][SPACE], 
               Boole B_symmetry_flag, 
               int   natom, 
               int   type[MAX_ATOMS] );
#endif
