#ifndef PRINTENERGIES
#define PRINTENERGIES
#include <stdio.h>
#include "autocomm.h"
void printEnergies( FloatOrDouble einter,
		    FloatOrDouble eintra,
		    FloatOrDouble torsFreeEnergy,
		    char  *prefixString, 
        int ligand_is_inhibitor
		    );
void printStateEnergies( FloatOrDouble einter,
			 FloatOrDouble eintra,
			 FloatOrDouble torsFreeEnergy,
			 char  *prefixString, 
			 int ligand_is_inhibitor
			 );
#endif
