
#ifndef EINTCALPRINT
#define EINTCALPRINT
#include "constants.h"
float  eintcalPrint( int   nonbondlist[MAX_NONBONDS][2], 
                float eint_table[NEINT][ATOM_MAPS][ATOM_MAPS], 
                float tcoord[MAX_ATOMS][SPACE], 
                int   atmtyp[MAX_ATOMS], 
                int   Nnb,
		Boole B_calcIntElec,
		float q1q2[MAX_NONBONDS]);
#endif
