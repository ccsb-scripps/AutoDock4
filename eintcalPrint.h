
#ifndef EINTCALPRINT
#define EINTCALPRINT
#include "constants.h"
FloatOrDouble  eintcalPrint( int   nonbondlist[MAX_NONBONDS][2], 
                FloatOrDouble eint_table[NEINT][ATOM_MAPS][ATOM_MAPS], 
                FloatOrDouble tcoord[MAX_ATOMS][SPACE], 
                int   atmtyp[MAX_ATOMS], 
                int   Nnb,
		Boole B_calcIntElec,
		FloatOrDouble q1q2[MAX_NONBONDS]);
#endif
