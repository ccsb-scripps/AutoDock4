
#ifndef EINTCAL
#define EINTCAL

#include "constants.h"

#ifndef EINTCALPRINT

FloatOrDouble  eintcal( int   nonbondlist[MAX_NONBONDS][4], 
                FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                FloatOrDouble tcoord[MAX_ATOMS][SPACE], 
                int   Nnb,
		Boole B_calcIntElec,
		FloatOrDouble q1q2[MAX_NONBONDS]);

#else	/*EINTCALPRINT*/

FloatOrDouble  eintcalPrint( int   nonbondlist[MAX_NONBONDS][4],
                FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                int   Nnb,
                Boole B_calcIntElec,
                FloatOrDouble q1q2[MAX_NONBONDS]);

#endif	/*EINTCALPRINT*/

#endif	/*!EINTCAL*/
